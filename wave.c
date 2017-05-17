#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <sys/time.h>

/* #define WRITE_TO_FILE */
//#define VERIFY

double global_error;
double total_time;

/** Define to enable debug mode */
#define DEBUG 0 /* 1 */

/** Debug output macro. Only active when DEBUG is non-0 */
#define dprintf(...)                    \
if (DEBUG)                              \
fprintf(stderr, __VA_ARGS__)

#define M_PI 3.141593
#define M_PI_2 1.570796
#define TRUE 1
#define FALSE 0

typedef struct {
    int cart_rank;/* reorder, row, column; */
    int px;
    int py;
    int x_coord;
    int y_coord;
    MPI_Comm cart_comm;
    
} info;

double timer();
double initialize(double x, double y, double t);
void save_solution(double *u, int Ny, int Nx, int n);
void makeGrid(info * grid);

int main(int argc, char *argv[])
{

    info grid;
   
    MPI_Status status;
    
    int Nx,Ny,Nt;
    double dt, dx, lambda_sq;
    double *u;
    double *u_old;
    double *u_new;
    double begin,end;
    
    Nx=128;
    if(argc>1)
    Nx=atoi(argv[1]);
    Ny=Nx;
    Nt=Nx;
    dx=1.0/(Nx-1);
    dt=0.50*dx;
    lambda_sq = (dt/dx)*(dt/dx);
    
    
    MPI_Init(&argc, &argv);
    
    makeGrid(&grid);
    
    
    int x_coord = grid.x_coord; int y_coord = grid.y_coord;
    
    int local_Nx = floor(Nx/grid.px);
    int local_Ny = floor(Ny/grid.py);
    
    int evenPx = 0; // check if dimensions are even
    int evenPy = 0;
    
    if (x_coord < Nx % grid.px) {
        local_Nx += 1;
        evenPx = 1;
    }
    if (y_coord < Ny % grid.py) {
        local_Ny += 1;
        evenPy = 1;
    }
    
    int xBlock = local_Nx + 2;
    int yBlock = local_Ny + 2;
    
    
    u = malloc(xBlock*yBlock*sizeof(double)); // add 1 to get extra row and column (halo points)
    u_old = malloc(xBlock*yBlock*sizeof(double));
    u_new = malloc(xBlock*yBlock*sizeof(double));
    
    
    /* Setup IC */
    memset(u,0,xBlock*yBlock*sizeof(double));       // BC
    memset(u_old,0,xBlock*yBlock*sizeof(double));   // BC
    memset(u_new,0,xBlock*yBlock*sizeof(double));   // BC
    
    
    
    // int MPI_Type_vector(int count, int blocklength, int stride, MPI_Datatype oldtype, MPI_Datatype *newtype)
    MPI_Datatype SEND_COLUMN;
    MPI_Type_vector(local_Ny, 1, xBlock, MPI_DOUBLE, &SEND_COLUMN);
    MPI_Type_commit(&SEND_COLUMN);
    
    
    MPI_Datatype SEND_ROW;
    MPI_Type_vector(local_Nx, 1, 1, MPI_DOUBLE, &SEND_ROW);
    MPI_Type_commit(&SEND_ROW);
    
    
    int ll_x = 1;
    int ll_y = 1;
    
    int ul_x = local_Nx + 1;
    int ul_y = local_Ny + 1;
    
    
    if (x_coord == 0)
    ll_x = 2;
    if (x_coord == grid.px - 1)
    ul_x = local_Nx;
    
    
    if (y_coord == 0)
    ll_y = 2;
    if (y_coord == grid.py - 1)
    ul_y = local_Ny;
    
    int x,y;
    for(int i = ll_y; i < ul_y; ++i) {
        for(int j = ll_x; j < ul_x; ++j) {
            
            
            if (evenPx == 1) {
                x = ((j-1) + local_Nx * x_coord) * dx;
            } else {
                x = ((j-1) + (Nx - (grid.px - x_coord) * local_Nx)) * dx;
            }
            
            if (evenPy == 1) {
                y = ((i-1) + local_Ny * y_coord) * dx;
            } else {
                y = ((i-1) + (Ny - (grid.py - y_coord) * local_Ny)) * dx;
            }
            u[y*xBlock + j] = initialize(x,y,0);
            u_new[y*xBlock + j] = initialize(x,y,dt);
        }
    }
    
#ifdef WRITE_TO_FILE
    save_solution(u_new,Ny,Nx,1); // write from one processor only
#endif
#ifdef VERIFY
    double max_error=0.0;
#endif
    
    /* Integrate */
    
    
    begin=timer();
    for(int n=2; n<Nt; ++n) {
        
        /* Swap ptrs */
        double *tmp = u_old;
        u_old = u;
        u = u_new;
        u_new = tmp;
        
        
           /* Communicate between processors */
        
        
        /* if not at left boundary, send left and recieve from left */
        if (x_coord > 0) {
            MPI_Send(&u[xBlock+1], 1, SEND_COLUMN, grid.cart_rank-1, 0, grid.cart_comm);
            MPI_Recv(&u[xBlock], 1, SEND_COLUMN, grid.cart_rank-1, 0, grid.cart_comm, &status);
        }
        
        /* if not at right boundary, send right and recieve from right */
        if (x_coord < grid.px - 1) {
            
            MPI_Recv(&u[xBlock + local_Nx + 1], 1, SEND_COLUMN, grid.cart_rank+1, 0, grid.cart_comm, &status);
            MPI_Send(&u[xBlock + local_Nx], 1, SEND_COLUMN, grid.cart_rank+1, 0, grid.cart_comm);
            
        }
        
        
        /* if not at botom, send down and recieve from bottom */
        if (y_coord < grid.py - 1) {
            
            MPI_Send(&u[xBlock * local_Ny +1], 1, SEND_ROW, grid.cart_rank + grid.px , 0, grid.cart_comm);
            MPI_Recv(&u[(local_Ny + 1) * xBlock + 1], 1, SEND_ROW,grid.cart_rank+grid.px, 0, grid.cart_comm, &status);
        }
        
        
        /* if not at top, send upwards and recieve from upper */
        if (y_coord > 0) {
            
            MPI_Recv(&u[1], 1, SEND_ROW, grid.cart_rank-grid.px, 0, grid.cart_comm, &status);
            MPI_Send(&u[xBlock + 1], 1, SEND_ROW, grid.cart_rank-grid.px, 0, grid.cart_comm);
        }
        
        
        /* Apply stencil */
        for(int i = ll_y + 1; i < ul_y - 1; ++i) {
            for(int j = ll_x + 1; j < ul_x - 1; ++j) {
                
                u_new[i*xBlock+j] = 2*u[i*xBlock+j]-u_old[i*xBlock+j]+lambda_sq*
                (u[(i+1)*xBlock+j] + u[(i-1)*xBlock+j] + u[i*xBlock+j+1] + u[i*xBlock+j-1] -4*u[i*xBlock+j]);
            }
        }
        
        
#ifdef VERIFY
        double error=0.0;
        
        
        for(int i = 1; i < local_Ny + 1; ++i) {
            for(int j = 1; j < local_Nx + 1; ++j) {
                
                if (evenPx == 1) {
                    x = ((j-1) + local_Nx * x_coord) * dx;
                } else {
                    x = ((j-1) + (Nx - (grid.px - x_coord) * local_Nx)) * dx;
                }
                
                
                if (evenPy == 1) {
                    y = ((i-1) + local_Ny * y_coord) * dx;
                } else {
                    y = ((i-1) + (Ny - (grid.py - y_coord) * local_Ny)) * dx;
                }
                
                
                double e = fabs(u_new[i*xBlock+j]-initialize(x,y,n*dt));
                if(e>error)
                    error = e;
            }
        }
        if(error > max_error)
            max_error=error;

#endif
        
#ifdef WRITE_TO_FILE
        save_solution(u_new,Ny,Nx,n);  // write from one processor only
#endif
        
    }
#ifdef VERIFY
    MPI_Reduce(&max_error, &global_error, 1, MPI_DOUBLE, MPI_MAX, 0, grid.cart_comm);
    
    if (grid.cart_rank == 0)
    printf("Maximum error: %g\n",global_error);
#endif
    end=timer();
    
    
    double local_time = end-begin;
    MPI_Reduce(&local_time, &total_time, 1, MPI_DOUBLE, MPI_MAX, 0, grid.cart_comm);
    
    
    if (grid.cart_rank == 0)
    
    printf("Time elapsed: %g s\n",total_time);
    
    free(u);
    free(u_old);
    free(u_new);
    MPI_Type_free(&SEND_ROW);
    MPI_Type_free(&SEND_COLUMN);
    
    
    MPI_Finalize();
    
    
    return 0;
}

double timer()
{
    struct timeval tv;
    gettimeofday(&tv, NULL);
    double seconds = tv.tv_sec + (double)tv.tv_usec / 1000000;
    return seconds;
}

double initialize(double x, double y, double t)
{
    double value = 0;
#ifdef VERIFY
    /* standing wave */
    value=sin(3*M_PI*x)*sin(4*M_PI*y)*cos(5*M_PI*t);
#else
    /* squared-cosine hump */
    const double width=0.1;
    
    double centerx = 0.25;
    double centery = 0.5;
    
    double dist = sqrt((x-centerx)*(x-centerx) +
                       (y-centery)*(y-centery));
    if(dist < width) {
        double cs = cos(M_PI_2*dist/width);
        value = cs*cs;
    }
#endif
    return value;
}

void save_solution(double *u, int Ny, int Nx, int n)
// TODO: enable saving solution when > 1 processor are used
{
    char fname[50];
    sprintf(fname,"solution-%d.dat",n);
    FILE *fp = fopen(fname,"w");
    
    fprintf(fp,"%d %d\n",Nx,Ny);
    
    for(int j = 0; j < Ny; ++j) {
        for(int k = 0; k < Nx; ++k) {
            fprintf(fp,"%e\n",u[j*Nx+k]);
        }
    }
    
    fclose(fp);
}
void makeGrid(info * grid) {
    int reorder = TRUE;
    int dims[2], period[2], coordinates[2];
    int nprocs;
    
    
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs); /* get current process id*/
    
    dims[0] = dims[1] = 0; // preventing bug making 1D typology
    
    MPI_Dims_create(nprocs, 2, dims);
    
    
    period[0] = period[1] = FALSE;
    
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, period, reorder, &(grid->cart_comm));
    MPI_Comm_rank(grid->cart_comm, &(grid->cart_rank));
    MPI_Cart_coords(grid -> cart_comm, grid -> cart_rank, 2, coordinates);
    
    grid -> x_coord = coordinates[1];
    grid -> y_coord = coordinates[0];
    
    grid -> py = dims[0];
    grid -> px = dims[1];
    
    
    
}



