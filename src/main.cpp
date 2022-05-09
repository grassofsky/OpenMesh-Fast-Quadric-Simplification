#include "fast_quadric_mesh_simplification.h"

#include <sstream>
#include <iostream>
#include <chrono>

void PrintMeshInfo(MyMesh& mesh)
{
    std::cout << "vertices: " << mesh.n_vertices() << std::endl;
    std::cout << "faces: " << mesh.n_faces() << std::endl;
}

void ShowHelp(const char* argv[])
{
    const char* cExe = argv[0];
    std::stringstream ssHelper;
    ssHelper << "Usage: " << cExe << " <input> <output> <vertex num of output>\n" << "\n"
        << "\tinput: the name of existing mesh file, supported file type is stl, off, obj, ply\n"
        << "\toutput: the name of decimated mesh file\n"
        << "\tvertex num of output: e.g. 1000\n";
    std::cout << ssHelper.str() << std::endl;
}

int main(int argc, const char *argv[])
{
    if (argc < 3)
    {
        ShowHelp(argv);
        return 0;
    }

    std::string strInput = argv[1];
    std::string strOutput = argv[2];
    int iVertexNum = atoi(argv[3]);

    MyMesh mesh;
    OpenMesh::IO::read_mesh(mesh, strInput);

    std::cout << "Before simplify: " << std::endl;
    PrintMeshInfo(mesh);

    std::chrono::system_clock::time_point tStart = std::chrono::system_clock::now();
    Simplification simplify(mesh);
    simplify.SimplifyVertexTo(iVertexNum);
    std::chrono::system_clock::time_point tEnd = std::chrono::system_clock::now();

    std::cout << "After simplify: " << std::endl;
    PrintMeshInfo(mesh);
    std::cout << "It took time: " << std::chrono::duration_cast<std::chrono::milliseconds>(tEnd - tStart).count() << " ms";

    OpenMesh::IO::write_mesh(mesh, strOutput);

    return 0;
}

