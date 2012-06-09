#include <stdio.h>
#include <string.h>
#include <mpi.h>
#include <math.h>

const int MAX_STRING=100;
double Trap (double, double, int, double);
double f(double);
void Get_input(int my_rank, int comm_sz, double* a_p, double*b_p, int* n_p);


int main(void)
{
	int		comm_sz;
	int		my_rank;
	int		n, local_n, i;
	double a, b, h, local_a, local_b;
	double local_int, total_int;
	int source;
	double local_traj[10];
	double global_traj[10][10];
	
	MPI_Init(NULL, NULL);
	MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

	//Get_input(my_rank, comm_sz, &a, &b, &n);
	
	/*h = (b-a)/n;
	local_n=n/comm_sz; 
	
	local_a = a + my_rank*local_n*h;
	local_b = local_a + local_n*h;
	local_int = Trap(local_a, local_b, local_n, h);
	*/
	
	if (my_rank != 0)
	{
		for (i = 0;i < 10; i++)
		{
			*(local_traj+i)=4.0;
			printf("%f \n",*(local_traj+i));
		}
		MPI_Send(&local_traj, 10, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
	}
	else 
	{	
		total_int = local_int; 
		for (source=1; source < comm_sz; source++) 
		{
			MPI_Recv(&local_traj, 10, MPI_DOUBLE, source, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			for (i = 0; i<10; i++)
			{
				//global_traj[source][i]=*(local_traj+i);
				//global_traj[source][i]=*(local_traj+sizeof(double));
				global_traj[source][i]=local_traj[i];	
				
				printf("%f \n",*(local_traj+i));

			}
		}
	}

	if (my_rank == 0)
	{
		for (source = 1; source < comm_sz; source++)
		{
			for (i = 0; i < 10; i++)
			{
			//	printf("With n = %d trapezoids, our estimation\n",n);
				printf("globalTraj[%d][%d]=%f \n",source, i, global_traj[source][i]);
	
			}
		}
	}
	MPI_Finalize();

	return 0;
}

double Trap( double a, double b, int n, double h)
{
	double estimate, x;
	int i;
	
	estimate = (f(a) + f(b))/2;
	for (i=1;i<=n;i++)
	{
		x = a+i*h;
		estimate+=f(x);
	}
	
	estimate=estimate*h;

}

void Get_input(int my_rank, int comm_sz, double* a_p, double*b_p, int* n_p)
{
	int dest;
	if (my_rank == 0)
	{
		printf("Enter a, b, and n\n");
		scanf("%lf %lf %d", a_p, b_p, n_p);
		for (dest = 1; dest < comm_sz; dest++)
		{
			MPI_Send(a_p, 1, MPI_DOUBLE, dest, 0, MPI_COMM_WORLD);
			MPI_Send(b_p, 1, MPI_DOUBLE, dest, 0, MPI_COMM_WORLD);
			MPI_Send(n_p, 1, MPI_INT, dest, 0, MPI_COMM_WORLD);
		}
	}
	else
	{
		MPI_Recv(a_p, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		MPI_Recv(b_p, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		MPI_Recv(n_p, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}
	return;
}



double f(double x)
{
	return log(x) ;
}