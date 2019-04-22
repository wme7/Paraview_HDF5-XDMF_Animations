#include <iostream>
#include <string>
#include <vector>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "hdf5.h"

// Compile in Linux as:
// g++ -std=c++11 -I /path/to/hdf5/include -L /path/to/hdf5/lib -lhdf5 test_HDF5_3DSMesh.cpp -o test3d.run
// In MacOS with homebrew HDF5 as:
// g++ -std=c++11 -I /usr/local/Cellar/hdf5/1.10.4/include -L /usr/local/Cellar/hdf5/1.10.4/lib -lhdf5 test_HDF5_3DSMesh.cpp -o test3d.run

void write_hdf5_geometry(
    const float *x,
    const float *y,
    const float *z,
    const uint nx,
    const uint ny,
    const uint nz)
{
    // Create file
    hid_t file_id;
    file_id = H5Fcreate("geometry.h5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    
    // Write the data file.
    hid_t dataset_id;
    hid_t dataspace_id;
    hsize_t Dims;
    hsize_t dims[3] = { nz, ny, nx };
    herr_t status;

    // Write number of cells in x-direction
    Dims = 1;
    dataspace_id = H5Screate_simple(1, &Dims, NULL);
    dataset_id = H5Dcreate2(file_id, "/NX_nodes", H5T_NATIVE_FLOAT, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &nx);
    status = H5Dclose(dataset_id);
    status = H5Sclose(dataspace_id);

    //  Write number of cells in y-direction
    Dims = 1;
    dataspace_id = H5Screate_simple(1, &Dims, NULL);
    dataset_id = H5Dcreate2(file_id, "/NY_nodes", H5T_NATIVE_FLOAT, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &ny);
    status = H5Dclose(dataset_id);
    status = H5Sclose(dataspace_id);

    //  Write number of cells in z-direction
    Dims = 1;
    dataspace_id = H5Screate_simple(1, &Dims, NULL);
    dataset_id = H5Dcreate2(file_id, "/NZ_nodes", H5T_NATIVE_FLOAT, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &nz);
    status = H5Dclose(dataset_id);
    status = H5Sclose(dataspace_id);

    // Write x coordinates
    dataspace_id = H5Screate_simple(1, &dims[2], NULL);
    dataset_id = H5Dcreate2(file_id, "/X_nodes", H5T_NATIVE_FLOAT, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, x);
    status = H5Dclose(dataset_id);
    status = H5Sclose(dataspace_id);

    // Write y coordinates
    dataspace_id = H5Screate_simple(1, &dims[1], NULL);
    dataset_id = H5Dcreate2(file_id, "/Y_nodes", H5T_NATIVE_FLOAT, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, y);
    status = H5Dclose(dataset_id);
    status = H5Sclose(dataspace_id);

    // Write z coordinates
    dataspace_id = H5Screate_simple(1, &dims[0], NULL);
    dataset_id = H5Dcreate2(file_id, "/Z_nodes", H5T_NATIVE_FLOAT, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, z);
    status = H5Dclose(dataset_id);
    status = H5Sclose(dataspace_id);

    // Close file
    status = H5Fclose(file_id);
}

void write_hdf5_field(
    const float *pressure,
    const uint nx,
    const uint ny,
    const uint nz,
    const uint it)
{
    std::string name = "fields_" + std::to_string(it) + ".h5";

    // Create file
    hid_t file_id;
    file_id = H5Fcreate(name.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    
    // Write the data file.
    hid_t dataset_id;
    hid_t dataspace_id;
    herr_t status;

    // Write the node centered data.
    hsize_t dims[3] = { nz, ny, nx };
    dataspace_id = H5Screate_simple(3, dims, NULL);
    dataset_id = H5Dcreate2(file_id, "/Pressure", H5T_NATIVE_FLOAT, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, pressure);
    status = H5Dclose(dataset_id);
    status = H5Sclose(dataspace_id);

    // Close file
    status = H5Fclose(file_id);
}
 
void write_xdmf_field(
    const uint nx,
    const uint ny,
    const uint nz,
    const uint it,
    const double t)
{
    FILE *xmf = 0;

    std::string name = "fields_" + std::to_string(it) + ".xmf";
 
    /*
     * Open the file and write the XML description of the mesh..
     */
    xmf = fopen(name.c_str(), "w");
    fprintf(xmf, "<?xml version=\"1.0\" ?>\n");
    fprintf(xmf, "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n");
    fprintf(xmf, "<Xdmf Version=\"2.0\">\n");
    fprintf(xmf, " <Domain>\n");
    fprintf(xmf, "   <Grid Name=\"mesh1\" GridType=\"Uniform\">\n");
    fprintf(xmf, "     <Topology TopologyType=\"3DRECTMesh\" NumberOfElements=\"%d %d %d\"/>\n", nz, ny, nx);
    fprintf(xmf, "     <Geometry GeometryType=\"VXVYVZ\">\n");
    fprintf(xmf, "       <DataItem Name=\"xcoords\" Dimensions=\"%d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n", nx);
    fprintf(xmf, "        geometry.h5:/X_nodes\n");
    fprintf(xmf, "       </DataItem>\n");
    fprintf(xmf, "       <DataItem Name=\"ycoords\" Dimensions=\"%d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n", ny);
    fprintf(xmf, "        geometry.h5:/Y_nodes\n");
    fprintf(xmf, "       </DataItem>\n");
    fprintf(xmf, "       <DataItem Name=\"zcoords\" Dimensions=\"%d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n", nz);
    fprintf(xmf, "        geometry.h5:/Z_nodes\n");
    fprintf(xmf, "       </DataItem>\n");
    fprintf(xmf, "     </Geometry>\n");
    fprintf(xmf, "     <Time TimeType=\"Single\" Value=\"%g\"/>\n", t);
    fprintf(xmf, "     <Attribute Name=\"Pressure\" AttributeType=\"Scalar\" Center=\"Node\">\n");
    fprintf(xmf, "       <DataItem Dimensions=\"%d %d %d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n", nz, ny, nx);
    fprintf(xmf, "        fields_%d.h5:/Pressure\n", it);
    fprintf(xmf, "       </DataItem>\n");
    fprintf(xmf, "     </Attribute>\n");
    fprintf(xmf, "   </Grid>\n");
    fprintf(xmf, " </Domain>\n");
    fprintf(xmf, "</Xdmf>\n");
    fclose(xmf);
}

void write_xdmf_temporal(
    const std::vector<size_t>& iterations, 
    const std::string name)
{
    FILE *xmf = 0;

    std::string name_xmf;
    name_xmf = name + ".xmf";

    /*
     * Open the file and write the XML description of the mesh.
     */
    xmf = fopen(name_xmf.c_str(), "w");
    fprintf(xmf, "<?xml version=\"1.0\" ?>\n");
    fprintf(xmf, "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n");
    fprintf(xmf, "<Xdmf xmlns:xi=\"http://www.w3.org/2001/XInclude\" Version=\"2.0\">\n");
    fprintf(xmf, "  <Domain>\n");
    fprintf(xmf, "   <Grid GridType=\"Collection\" CollectionType=\"Temporal\">\n");
    for(auto it : iterations)
    {
        fprintf(xmf, "    <xi:include href=\"fields_%d.xmf\"  xpointer=\"xpointer(//Xdmf/Domain/Grid)\"/>\n", int(it));
    }
    fprintf(xmf, "   </Grid>\n");
    fprintf(xmf, "  </Domain>\n");
    fprintf(xmf, "</Xdmf>\n");
    fclose(xmf);
}
 
int main(int argc, char *argv[])
{
    // Parameters
    double dt=0.1;
    double t=0;
    std::string name = "numSolution";

    // HDF5 convention is fastest is last, this is {nz, ny, nx}.
    hsize_t nx = 256;
    hsize_t ny = 192;
    hsize_t nz = 128;
   
    // Arrays for grid spacing along each axis.
    float ord = 0.0f;
    std::vector<float> xcoords;
    for (int i = 0; i < nx; i++) {
        xcoords.push_back(ord);
        ord += i * 0.1f;
    }

    ord = 0.0f;
    std::vector<float> ycoords;
    for (int i = 0; i < ny; i++) {
        ycoords.push_back(ord);
        ord += i * 0.1f;
    }

    ord = 0.0f;
    std::vector<float> zcoords;
    for (int i = 0; i < nz; i++) {
        zcoords.push_back(ord);
        ord += i * 0.1f;
    }
    // Save the geometry of the problem
    write_hdf5_geometry(&xcoords.front(),&ycoords.front(),&zcoords.front(),nx,ny,nz);

    // Iteration loop
    std::vector<size_t> iterations;
    for (size_t it = 0; it < 10; it++)
    {
        // iteration time 
        t = dt*it;

        // arbitrary scalar data for writing.
        std::vector<float> scalars;
        float total = 256 * 192;
        for (int k = 0; k < nz; k++) {
            for (int j = 0; j < ny; j++) {
                for (int i = 0; i < nx; i++) {
                    float x = float(i) / (atan(1.0) * nx);
                    float y = 2.0f * float(j) / ny;
                    scalars.push_back( sinf(x) * cosf(y) * (1.0-t) );
                }
            }
        }
        // Save fields information
        write_hdf5_field(&scalars.front(),nx,ny,nz,it);
        write_xdmf_field(nx,ny,nz,it,t);
        // accumulate iteration history
        iterations.push_back(it);
    }
    // Lastly create the iteration fields list
    write_xdmf_temporal(iterations,name);
 
    return 0;
}