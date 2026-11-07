#include "fvCFD.H"
#include "OFstream.H"
int main(int argc,char* argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    fileName out(runTime.path()/"constant"/"cell_geometry.csv");
    OFstream os(out);
    os.precision(15);
    os<<"x_m,y_m,z_m,volume_m3\n";
    forAll(mesh.C(),i)
    {
        const vector& p=mesh.C()[i];
        os<<p.x()<<","<<p.y()<<","<<p.z()<<","<<mesh.V()[i]<<"\n";
    }
    Info<<"Wrote "<<mesh.nCells()<<" cell geometries to "<<out<<endl;
    return 0;
}
