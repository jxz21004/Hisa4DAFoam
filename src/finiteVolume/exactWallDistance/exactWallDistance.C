// Nearest distance to the triangulated global wall surface.
// Opt-in runtime plugin; no change to the default OpenFOAM wall-distance method.
#include "patchDistMethod.H"
#include "fvMesh.H"
#include "volFields.H"
#include "triSurface.H"
#include "triSurfaceSearch.H"
#include "Pstream.H"
#include "addToRunTimeSelectionTable.H"
#include <vector>
#include <cmath>

namespace Foam { namespace patchDistMethods {
class exactWallDistance : public patchDistMethod
{
    List<labelList> associations_;
    std::vector<double> wallSnapshot_;
    List<std::vector<double>> querySnapshots_;
    static double passive(const scalar& x)
    {
#if defined(CODI_ADR) || defined(CODI_ADF)
        return x.getValue();
#else
        return x;
#endif
    }
    static std::vector<double> snapshot(const pointField& points)
    {
        std::vector<double> values; values.reserve(3*points.size());
        for (const point& p : points) for (direction j=0;j<3;++j) values.push_back(passive(p[j]));
        return values;
    }
    bool calculate(volScalarField& y, volVectorField* normals)
    {
        DynamicList<point> local;
        for (const label patchi : patchIDs_.sortedToc())
        {
            const polyPatch& patch = mesh_.boundaryMesh()[patchi];
            forAll(patch, facei)
            {
                const face& f=patch[facei];
                // Match face::nearestPoint geometry, including warped polygons.
                if (f.size()==3)
                {
                    for (const label vertex : f) local.append(mesh_.points()[vertex]);
                }
                else
                {
                    const point centre=f.centre(mesh_.points());
                    forAll(f,i)
                    {
                        local.append(mesh_.points()[f[i]]);
                        local.append(mesh_.points()[f[(i+1)%f.size()]]);
                        local.append(centre);
                    }
                }
            }
        }
        // Communicate typed scalar buffers: nested Field serialization can
        // discard reverse-mode dependencies on remote wall coordinates.
        labelList counts(Pstream::nProcs(),0);
        counts[Pstream::myProcNo()]=local.size();
        Pstream::gatherList(counts); Pstream::broadcastList(counts);
        label count=0, start=0;
        forAll(counts,rank)
        {
            if (rank<Pstream::myProcNo()) start+=counts[rank];
            count+=counts[rank];
        }
        if (!count || count%3)
            FatalErrorInFunction << "Empty or invalid global wall triangulation" << exit(FatalError);
        scalarField packed(3*count,Zero);
        forAll(local,i) for(direction j=0;j<3;++j) packed[3*(start+i)+j]=local[i][j];
        if (Pstream::parRun()) UPstream::mpiAllReduce
        (
            packed.data(),packed.size(),UPstream::opCodes::op_sum,UPstream::worldComm
        );
        pointField points(count);
        forAll(points,i) for(direction j=0;j<3;++j) points[i][j]=packed[3*i+j];
        List<labelledTri> triangles(count/3);
        forAll(triangles,i) triangles[i]=labelledTri(3*i,3*i+1,3*i+2,0);
        const auto currentWall=snapshot(points);
        const bool sameWall=currentWall==wallSnapshot_;
        const label nSets=1+mesh_.boundary().size();
        if (associations_.size()!=nSets)
        {
            associations_.setSize(nSets);
            querySnapshots_.setSize(nSets);
        }
        autoPtr<triSurface> surface;
        autoPtr<triSurfaceSearch> search;
        label seti=0, reused=0, searched=0;
        auto sample = [&](const pointField& centres, scalarField& distance, vectorField* n)
        {
            labelList& ids=associations_[seti];
            const auto currentQuery=snapshot(centres);
            const bool reuse=sameWall && ids.size()==centres.size()
                && querySnapshots_[seti]==currentQuery;
            if (!reuse)
            {
                if (!search.valid())
                {
                    surface.reset(new triSurface(triangles, geometricSurfacePatchList(), points));
                    search.reset(new triSurfaceSearch(surface()));
                }
                scalarField bounds(centres.size(),GREAT);
                // Old associations supply an upper bound, not a frozen nearest face.
                if (ids.size()==centres.size()) forAll(centres,i)
                {
                    const label id=ids[i];
                    if (id>=0 && id<triangles.size())
                    {
                        const triPointRef tri(points[3*id],points[3*id+1],points[3*id+2]);
                        const scalar d=tri.nearestPoint(centres[i]).distance();
                        bounds[i]=sqr(d)+max(sqr(d)*scalar(1e-10),scalar(SMALL));
                    }
                }
                List<pointIndexHit> hits;
                search().findNearest(centres,bounds,hits);
                ids.setSize(centres.size());
                forAll(centres,i)
                {
                    if (!hits[i].hit()) FatalErrorInFunction
                        << "No wall hit for " << centres[i] << exit(FatalError);
                    ids[i]=hits[i].index();
                }
                searched+=centres.size();
            }
            else reused+=centres.size();
            // Always re-evaluate projection with current active coordinates.
            // Never cache AD scalars or numerical distances across recordings.
            forAll(centres,i)
            {
                const label id=ids[i];
                const triPointRef tri(points[3*id],points[3*id+1],points[3*id+2]);
                const vector delta=tri.nearestPoint(centres[i]).rawPoint()-centres[i];
                distance[i]=mag(delta);
                if (n) (*n)[i]=delta/max(distance[i],VSMALL);
            }
            querySnapshots_[seti]=currentQuery;
            ++seti;
        };
        sample(mesh_.C().primitiveField(), y.primitiveFieldRef(),
            normals ? &normals->primitiveFieldRef() : nullptr);
        forAll(mesh_.boundary(),patchi)
        {
            auto& yp=y.boundaryFieldRef()[patchi];
            if (patchIDs_.found(patchi))
            {
                // Match meshWave boundary convention: SA evaluates 0/d^2 on walls.
                yp == scalar(SMALL);
                ++seti;
                if (normals) normals->boundaryFieldRef()[patchi] == mesh_.boundary()[patchi].nf();
            }
            else
            {
                vectorField* np=normals ? &normals->boundaryFieldRef()[patchi] : nullptr;
                sample(mesh_.boundary()[patchi].Cf(), yp, np);
                forAll(yp,i) yp[i] += SMALL;
            }
        }
        // Keep OpenFOAM patch coupling/BC semantics.
        y.correctBoundaryConditions();
        if (normals) normals->correctBoundaryConditions();
        Info<< "exactWallDistance: " << triangles.size()
            << " global triangles; cell distance min/max "
            << gMin(y.primitiveField()) << ' ' << gMax(y.primitiveField()) << endl;
        wallSnapshot_=currentWall;
        Info<< "exactWallDistance associations: reused=" << returnReduce(reused,sumOp<label>())
            << " searched=" << returnReduce(searched,sumOp<label>()) << endl;
        return false;
    }
public:
    TypeName("exactWallDistance");
    exactWallDistance(const dictionary&,const fvMesh& mesh,const labelHashSet& ids)
    : patchDistMethod(mesh,ids) {}
    void updateMesh(const mapPolyMesh&) override
    {
        associations_.clear(); querySnapshots_.clear(); wallSnapshot_.clear();
    }
    bool correct(volScalarField& y) override { return calculate(y,nullptr); }
    bool correct(volScalarField& y,volVectorField& n) override { return calculate(y,&n); }
};
defineTypeNameAndDebug(exactWallDistance,0);
addToRunTimeSelectionTable(patchDistMethod,exactWallDistance,dictionary);
}}
