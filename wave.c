#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <sys/time.h>

/* #define WRITE_TO_FILE */
#define VERIFY



/** Define to enable debug mode */
#define DEBUG 1 /* 1 */

/** Debug output macro. Only active when DEBUG is non-0 */
#define dprintf(...)                    \
if (DEBUG)                              \
fprintf(stderr, __VA_ARGS__)

//#define M_PI 3.141593
//#define M_PI_2 1.570796
#define TRUE 1
#define FALSE 0

typedef struct {
    int cart_rank, nprocs, reorder, row, column;
    int px;
    int py;
    MPI_Comm cart_comm, row_comm, col_comm;
    
} info;


double timer();
double initialize(double x, double y, double t);
void save_solution(double *u, int Ny, int Nx, int n);
void makeGrid(info * grid);

int main(int argc, char *argv[])
{
    
    MPI_Init(&argc, &argv);
    
    int global_rank;
    info grid;
    int nprocs;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &global_rank);
    
    
    grid.nprocs = nprocs;
    makeGrid(&grid);
    
    

    if (global_rank == 2) {
    
        dprintf("px %d \n", grid.px);
        dprintf("py %d \n", grid.py);
        
    }
    
    
    MPI_Finalize();
    
    
  int Nx,Ny,Nt, N;
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
  N = Nx;
  dx=1.0/(Nx-1);
  dt=0.50*dx;
  lambda_sq = (dt/dx)*(dt/dx);

  u = malloc(Nx*Ny*sizeof(double));
  u_old = malloc(Nx*Ny*sizeof(double));
  u_new = malloc(Nx*Ny*sizeof(double));

  /* Setup IC */

  memset(u,0,Nx*Ny*sizeof(double));
  memset(u_old,0,Nx*Ny*sizeof(double));
  memset(u_new,0,Nx*Ny*sizeof(double));

  for(int i = 1; i < (Ny-1); ++i) {
    for(int j = 1; j < (Nx-1); ++j) {
      double x = j*dx;
      double y = i*dx;

      /* u0 */
      u[i*Nx+j] = initialize(x,y,0);

      /* u1 */
      u_new[i*Nx+j] = initialize(x,y,dt);
    }
  }

#ifdef WRITE_TO_FILE
  save_solution(u_new,Ny,Nx,1);
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

    /* Apply stencil */
    for(int i = 1; i < (Ny-1); ++i) {
      for(int j = 1; j < (Nx-1); ++j) {

        u_new[i*Nx+j] = 2*u[i*Nx+j]-u_old[i*Nx+j]+lambda_sq*
          (u[(i+1)*Nx+j] + u[(i-1)*Nx+j] + u[i*Nx+j+1] + u[i*Nx+j-1] -4*u[i*Nx+j]);
      }
    }

#ifdef VERIFY
    double error=0.0;
    for(int i = 0; i < Ny; ++i) {
      for(int j = 0; j < Nx; ++j) {
        double e = fabs(u_new[i*Nx+j]-initialize(j*dx,i*dx,n*dt));
        if(e>error)
          error = e;
      }
    }
    if(error > max_error)
      max_error=error;
#endif

#ifdef WRITE_TO_FILE
    save_solution(u_new,Ny,Nx,n);
#endif

  }
  end=timer();

  printf("Time elapsed: %g s\n",(end-begin));

#ifdef VERIFY
  printf("Maximum error: %g\n",max_error);
#endif

  free(u);
  free(u_old);
  free(u_new);

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
    int reorder = FALSE;
    int dims[2], period[2], coordinates[2], row_col_coordinates[2];
    int global_rank;
    
    MPI_Comm_size(MPI_COMM_WORLD, &(grid->nprocs)); /* get current process id*/
    MPI_Comm_rank(MPI_COMM_WORLD, &global_rank); /* get current process id */
    
    /* int MPI_Dims_create(int nnodes, int ndims, int dims[]) */
    
    dims[0] = dims[1] = 0; // preventing bug making 1D typology
    
    MPI_Dims_create(grid -> nprocs, 2, dims);

    
//    if (global_rank == 1) {
//        
//        dprintf("px %d \n", dims[0]);
//        dprintf("py %d \n", dims[1]);
//    }
    
  //  dim[0] = dim[1] = sqrt((grid->nprocs)); /*dimension*/
    period[0] = period[1] = TRUE;
    
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, period, reorder, &(grid->cart_comm));
    MPI_Comm_rank(grid->cart_comm, &(grid->cart_rank));
    MPI_Cart_coords(grid->cart_comm, grid->cart_rank, 2, coordinates);
    
    grid -> row = coordinates[0];
    grid -> column = coordinates[1];
    //grid->q = sqrt(grid->nprocs);
    
    grid -> px = dims[0];
    grid -> py = dims[1];
   
   
    /*Create subcommunicators for row and column*/
    row_col_coordinates[0] = 0;
    row_col_coordinates[1] = 1;
    MPI_Cart_sub(grid->cart_comm, row_col_coordinates, &(grid->row_comm)); //row communicator
    row_col_coordinates[0] = 1;
    row_col_coordinates[1] = 0;
    MPI_Cart_sub(grid->cart_comm, row_col_coordinates, &(grid->col_comm)); //column communicator
    
    
}


