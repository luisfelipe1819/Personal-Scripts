/*DICCIONARIO
conexn.-Variable temporal que guarda el numero ingresado por el usuario cuando se preguntan las -conexiones- entre los nodos.
mRed[][].-Matriz en donde se guarda la red(-matriz de la Red-).
i.-Contador.
j.-Contador.
k.-Contador.
tnod.-Variable en donde se guarda el -total de nodos- de la red.
gNod.-Variable donde se guarda temporalmente el grado del nodo actual.
vClost[].- Vector donde se guardan temporalmente los nodos con los cuales tiene conexion nuestro nodo actual.
conexEntreNod.-Variable donde se guardan momentaneamente el numero de conexiones entre los nodos que tienen conexion con nuestro nodo actual. Es float para que el resultado de la operacion del coefifiente de clustering sean flotante.

*/

#include <stdio.h>

int main () { 

	short int conexn, mRed[117][117], i, j, tnod, k, gNod, vClost[116];

	float conexEntreNod;

	printf ("\n\n¿Cuantos nodos tiene la red? (Maximo 117 nodos)\n\n");
/*Se ingresa el total de nodos de la red*/
	scanf ("%d", &tnod);

	for (i = 0; i < tnod; i++)                                     /**/

		for (j = i; j < tnod; j++) {                               /**/
		
			if (j == i)                                             /**/

				mRed[i][j] = 0;                                   /**/

			else {                                               /**/
			
				printf ("\n\n¿Hay conexion entre el nodo %d y el nodo %d?(Si = 1, No = 0)\n\n", i + 1, j + 1);

				scanf ("%d", &conexn);                          /**/

				if ( conexn == 1 || conexn == 0 ) {              /**/
				
					mRed[i][j] = conexn;                        /*Se pregunta al usuario si hay conexion entre cada par de nodos sin repetirlos. El usuario debe de ingresar 0 y 1, si no lo hace salta un mensaje de error y le vuelve a preguntar la misma conexion.*/

					mRed[j][i] = conexn;                        /**/
				
				}                                                /**/

				else {                                           /**/
				
				printf ("\n\nIngrese caracter valido\n\n");      /**/

				j--;                                             /**/
				
				}                                                /**/
			
			}                                                    /**/
		
		}                                                        /**/

	printf ("La red :\n\t");               /**/

	for (i = 0; i < tnod; i++)             /**/

		printf ("[%d]\t", i + 1);          /**/

	for (i = 0; i < tnod; i++) {           /*Se imprime la red que ingreso el usuario*/

		printf ("\n\n[%d]\t", i + 1);      /**/

		for (j = 0; j < tnod; j++)          /**/

			printf ("%d\t", mRed[i][j]);    /**/

	}                                       /**/

	for (i = 0; i < tnod; i++) {
	
		printf ("\n\nEl nodo %d\t", i + 1);    /*Se va imprimiendo parte del mensaje que indica el nodo, su grado y su coeficiente*/

		gNod = 0;           /*Se inicializa en cero para evitar errores*/

		conexEntreNod = 0;  /*Se inicializa en cero para evitar errores*/

		for (j = 0; j < tnod; j++)  /**/
		
			if ( mRed[i][j])        /*Recorre todas las posibles conexiones del nodo y en caso de encontrarlas las va sumando*/

				vClost[gNod++] = j; /*Guarda en un vector todos los nodos con los que se van encontrando que tiene conexion el nodo actual*/
		
		printf (" Tiene un grado %d\t", gNod); /*Imprime el grado del nodo*/

		for ( j = 0; j < gNod; j++ )              /**/

			for ( k = j; k < gNod; k++ )          /*Se evaluan todas la conexiones entre los nodos que tuvo conexion nuestro nodo actual son repetirlas. Posteriormente las guarda*/
			
				if ( mRed[vClost[j]][vClost[k]] ) /**/
			
					conexEntreNod++;              /**/
		
		if (gNod == 1 || gNod == 0)

			printf (" \nY un coeficiente de clustering de 0.000000\t");

		else

		printf (" \nY un coeficiente de clustering de %.6f\t", (conexEntreNod * 2) / (gNod * (gNod - 1))); /*Se calcula el coeficiente de clustering con una presicion de 6 digitos decimales y se imprime*/
	
	}

	return 0;

}
