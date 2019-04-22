#include <iostream>
#include <string>
#include <vector>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "hdf5.h"

// Compile Linux as:
// g++ -std=c++11 -I /path/to/hdf5/include -L /path/to/hdf5/lib -lhdf5 test_HDF5_2DSMesh.cpp -o test2d.run
// In MacOS with homebrew HDF5 as:
// g++ -std=c++11 -I /usr/local/Cellar/hdf5/1.10.4/include -L /usr/local/Cellar/hdf5/1.10.4/lib -lhdf5 test_HDF5_2DSMesh.cpp -o test2d.run 

void write_hdf5_geometry(
    const float *x,
    const float *y,
    const uint nx,
    const uint ny)
{
    // Create file
    hid_t file_id;
    file_id = H5Fcreate("geometry.h5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    
    // Write the data file.
    hid_t dataset_id;
    hid_t dataspace_id;
    hsize_t dim;
    hsize_t coordDims[2];
    int inputbuffer;
    herr_t status;

    // Write number of cells in x-direction
    dim = 1;
    dataspace_id = H5Screate_simple(1, &dim, NULL);
    dataset_id = H5Dcreate2(file_id, "/NX_cells", H5T_NATIVE_FLOAT, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &nx);
    status = H5Dclose(dataset_id);
    status = H5Sclose(dataspace_id);

    //  Write number of cells in y-direction
    dim = 1;
    dataspace_id = H5Screate_simple(1, &dim, NULL);
    dataset_id = H5Dcreate2(file_id, "/NY_cells", H5T_NATIVE_FLOAT, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &ny);
    status = H5Dclose(dataset_id);
    status = H5Sclose(dataspace_id);

    // Write number of cells in x-direction
    dim = 1;
    inputbuffer = nx+1;
    dataspace_id = H5Screate_simple(1, &dim, NULL);
    dataset_id = H5Dcreate2(file_id, "/NX_nodes", H5T_NATIVE_FLOAT, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &inputbuffer);
    status = H5Dclose(dataset_id);
    status = H5Sclose(dataspace_id);

    //  Write number of cells in y-direction
    dim = 1;
    inputbuffer = ny+1;
    dataspace_id = H5Screate_simple(1, &dim, NULL);
    dataset_id = H5Dcreate2(file_id, "/NY_nodes", H5T_NATIVE_FLOAT, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &inputbuffer);
    status = H5Dclose(dataset_id);
    status = H5Sclose(dataspace_id);

    // Write x coordinates
    coordDims[0] = (ny + 1);
    coordDims[1] = (nx + 1);
    dataspace_id = H5Screate_simple(2, coordDims, NULL);
    dataset_id = H5Dcreate2(file_id, "/X_nodes", H5T_NATIVE_FLOAT, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, x);
    status = H5Dclose(dataset_id);
    status = H5Sclose(dataspace_id);

    // Write y coordinates
    coordDims[0] = (ny + 1);
    coordDims[1] = (nx + 1);
    dataspace_id = H5Screate_simple(2, coordDims, NULL);
    dataset_id = H5Dcreate2(file_id, "/Y_nodes", H5T_NATIVE_FLOAT, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, y);
    status = H5Dclose(dataset_id);
    status = H5Sclose(dataspace_id);

    // Close file
    status = H5Fclose(file_id);
}

void write_hdf5_field(
    const float *pressure,
    const float *velocityx,
    const uint nx,
    const uint ny,
    const uint it)
{
    std::string name = "fields_" + std::to_string(it) + ".h5";

    // Create file
    hid_t file_id;
    file_id = H5Fcreate(name.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    
    // Write the data file.
    hid_t dataset_id;
    hid_t dataspace_id;
    hsize_t dims[2];
    herr_t status;
   
    // Write the cell centered data.
    dims[0] = ny;
    dims[1] = nx;
    dataspace_id = H5Screate_simple(2, dims, NULL);
    dataset_id = H5Dcreate2(file_id, "/Pressure", H5T_NATIVE_FLOAT, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, pressure);
    status = H5Dclose(dataset_id);
    status = H5Sclose(dataspace_id);

    // Write the node centered data.
    dims[0] = ny + 1;
    dims[1] = nx + 1;
    dataspace_id = H5Screate_simple(2, dims, NULL);
    dataset_id = H5Dcreate2(file_id, "/VelocityX", H5T_NATIVE_FLOAT, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, velocityx);
    status = H5Dclose(dataset_id);
    status = H5Sclose(dataspace_id);

    // Close file
    status = H5Fclose(file_id);
}
 
void write_xdmf_field(
    const uint nx,
    const uint ny,
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
    fprintf(xmf, "     <Topology TopologyType=\"2DSMesh\" NumberOfElements=\"%d %d\"/>\n", ny+1, nx+1);
    fprintf(xmf, "     <Geometry GeometryType=\"X_Y\">\n");
    fprintf(xmf, "       <DataItem Dimensions=\"%d %d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n", (ny+1), (nx+1));
    fprintf(xmf, "        geometry.h5:/X_nodes\n");
    fprintf(xmf, "       </DataItem>\n");
    fprintf(xmf, "       <DataItem Dimensions=\"%d %d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n", (ny+1), (nx+1));
    fprintf(xmf, "        geometry.h5:/Y_nodes\n");
    fprintf(xmf, "       </DataItem>\n");
    fprintf(xmf, "     </Geometry>\n");
    fprintf(xmf, "     <Time TimeType=\"Single\" Value=\"%g\"/>\n", t);
    fprintf(xmf, "     <Attribute Name=\"Pressure\" AttributeType=\"Scalar\" Center=\"Cell\">\n");
    fprintf(xmf, "       <DataItem Dimensions=\"%d %d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n", ny, nx);
    fprintf(xmf, "        fields_%d.h5:/Pressure\n", it);
    fprintf(xmf, "       </DataItem>\n");
    fprintf(xmf, "     </Attribute>\n");
    fprintf(xmf, "     <Attribute Name=\"VelocityX\" AttributeType=\"Scalar\" Center=\"Node\">\n");
    fprintf(xmf, "       <DataItem Dimensions=\"%d %d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n", ny+1, nx+1);
    fprintf(xmf, "        fields_%d.h5:/VelocityX\n", it);
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
    // The number of cells in the X, Y dimensions
    int nx=30;
    int ny=20;
    double dt=0.1;
    double t=0;
    std::vector<size_t> iterations;
    std::string name = "numSolution";
    // #define M_PI 3.1415926535897932

    // Create the coordinate data.
    float *x = (float *) malloc((nx+1)*(ny+1) * sizeof(float));
    float *y = (float *) malloc((nx+1)*(ny+1) * sizeof(float));
    int ndx = 0;
    for (int j = 0; j < ny+1; j++)
    {
        float yt = j*1.0 / ny;
        float angle = yt * M_PI;
        for (int i = 0; i < nx+1; i++)
        {
            float xt = i*1.0 / nx;
            float R = (1.-xt)*2. + xt*5.;
 
            x[ndx] = R * cos(angle);
            y[ndx] = R * sin(angle);
            ndx++;
        }
    }
    // Save the geometry of the problem
    write_hdf5_geometry(x,y,nx,ny);

    // Create the scalar data
    float *pressure = (float *) malloc(nx*ny * sizeof(float));
    float *velocityx = (float *) malloc((nx+1)*(ny+1) * sizeof(float));
    // Iteration loop
    for (size_t it = 0; it < 20; it++)
    {
        // iteration time 
        t = dt*it;

        // Create the scalar data for pressure 
        for (int j = 0; j < ny; j++)
        {
            for (int i = 0; i < nx; i++)
            {
                int ndx = j * nx + i;
                pressure[ndx] = float(j+2*t);
            }
        }
        // Create the scalar data for velocityx
        for (int j = 0; j < ny+1; j++)
        {
            for (int i = 0; i < nx+1; i++)
            {
                int ndx = j * (nx+1) + i;
                velocityx[ndx] = float(i+10*t);
            }
        }
        // Save fields information
        write_hdf5_field(pressure,velocityx,nx,ny,it);
        write_xdmf_field(nx,ny,it,t);
        // accumulate iteration history
        iterations.push_back(it);
    }
    // Lastly create the iteration fields list
    write_xdmf_temporal(iterations,name);

    // Free the data.
    free(x);
    free(y);
    free(pressure);
    free(velocityx);
 
    return 0;
}