#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//#include <openblas/lapacke.h>
#include <mkl_lapacke.h>

double *generate_matrix(int size)
{
  int i;
  double *matrix = (double *) malloc(sizeof(double) * size * size);

  srand(1);

  for (i = 0; i < size * size; i++) {
    matrix[i] = rand() % 100;
  }

  return matrix;
}

int is_nearly_equal(double x, double y)
{
  const double epsilon = 1e-5 /* some small number */;
  return abs(x - y) <= epsilon * abs(x);
  // see Knuth section 4.2.2 pages 217-218
}

int check_result(double *bref, double *b, int size)
{
  int i;

  for(i = 0; i < size*size; i++) {
    if (!is_nearly_equal(bref[i], b[i]))
      return 0;
  }

  return 1;
}

int my_dgesv(int n, int nrhs, double *a, int lda, int *ipiv, double *b, int ldb)
{

        double identidad[n][n];
	double matrix_a[n][n];
	double matrix_b[n][n];
	double matrix_x[n][n];

	int i, j, k, l;

	//creamos una matriz identidad                                                     
	for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            if (i == j) {
                identidad[i][j] = 1;
                }
                else {
                        identidad[i][j]=0;
                } 
          }
        }
	
  //pasamos las matrices A y B al formato que usamos para la identidad                      
    k = 0;
    for (i = 0; i < n; i++) {
      k = i*n;
        for (j = 0; j < n; j++) {
	    matrix_a[i][j] = a[k];
            matrix_b[i][j] = b[k];
	    k++;
        }
    }	

    //recorremos los primeros elementos de cada fila para ponerlos a 0 hasta alcanzar el central, la diagonal pondremos a 1           
    for (i = 0; i < n; i++){
    	for (j = 0; j < n; j++){
            if (i>j){
	      double variable;
	      int indice;
             //buscamos la fila de la matriz que vale 1 en la misma columna
	      indice = j;
	       //calculamos el numero por el que se multiplicara toda la fila
                variable = matrix_a[i][j]/matrix_a[indice][j];
                double array_matriz[n], array_identidad[n];
                //restamos la fila multiplicada por el valor de la fila actual
                for (l = 0; l < n; l++){
                	array_matriz[l] = matrix_a[indice][l]*variable;
                	matrix_a[i][l] = matrix_a[i][l] - array_matriz[l];
                	array_identidad[l] = identidad[indice][l]*variable;
                	identidad[i][l] = identidad[i][l] - array_identidad[l];
				}
            }
            //ponemos los elementos centrales en valor 1 y convertimos toda la fila                          
            else if (i == j){
            	double variable;
            	variable = matrix_a[i][j];
            	for(k = 0; k < n; k++){
            		matrix_a[i][k] = matrix_a[i][k]/variable;
			identidad[i][k] = identidad[i][k]/variable;
       			}
	    	}
        }
    }

    //recorremos la matriz desde la esquina inferior derecha para poner a 0 todo a la derecha de la diagonal 
    for (i = (n-1); i >= 0; i--){
    	for (j = (n-1); j >= 0; j--){
            if (i<j){
	      double variable;
	      int indice, l;
              //buscamos la fila de la matriz que vale 1 en la misma columna
	      indice = j;
			   	
	     	//calculamos el numero por el que se multiplicara toda la fila
                variable = matrix_a[i][j]/matrix_a[indice][j];
                double array_matriz[n], array_identidad[n];
                //restamos la fila multiplicada por el valor de la fila actual
                for (l = (n-1); l >= 0; l--){
                	array_matriz[l] = matrix_a[indice][l]*variable;
                	matrix_a[i][l] = matrix_a[i][l] - array_matriz[l];
                	array_identidad[l] = identidad[indice][l]*variable;
                	identidad[i][l] = identidad[i][l] - array_identidad[l];
				}
            }
        }
    }
   
    //multiplicamos A^-1 * B para hallar X
    for (i = 0; i < n; i++){
    	for (j = 0; j < n; j++){
		double valor = 0;
		//suma de fila por columna
		for(k = 0; k < n; k++){
		  valor = valor + (identidad[i][k] * matrix_b[k][j]);
		}
		matrix_x[i][j] = valor;
	}
    }
   
    k = 0;
    for(i = 0; i < n; i++){
      k = i*n;
      for(j = 0; j < n; j++){
	      b[k] = matrix_x[i][j];
	      k++;
      }
    }
  
}

void main(int argc, char *argv[])
{
  int size = atoi(argv[1]);

  double *a, *aref;
  double *b, *bref;

  a = generate_matrix(size);
  aref = generate_matrix(size);
  b = generate_matrix(size);
  bref = generate_matrix(size);

  // Using LAPACK dgesv OpenBLAS implementation to solve the system
  int n = size, nrhs = size, lda = size, ldb = size, info;
  int *ipiv = (int *) malloc(sizeof(int) * size);

  clock_t tStart = clock();
  info = LAPACKE_dgesv(LAPACK_ROW_MAJOR, n, nrhs, aref, lda, ipiv, bref, ldb);
  printf("Time taken by OpenBLAS LAPACK: %.2fs\n", (double) (clock() - tStart) / CLOCKS_PER_SEC);

  int *ipiv2 = (int *) malloc(sizeof(int) * size);

  tStart = clock();
  my_dgesv(n, nrhs, a, lda, ipiv2, b, ldb);
  printf("Time taken by my implementation: %.2fs\n", (double) (clock() - tStart) / CLOCKS_PER_SEC);
    
  if (check_result(bref, b, size) == 1)
    printf("Result is ok!\n");
  else
    printf("Result is wrong!\n");
}
