/*
    Proyecto final de la materia Metodos numericos en la ingenieria
    Integrantes:
    Alban Aguilar Campos
    Noe Flores Sifuentes
    Salvador Garcia Martinez
    Francisco Navarrete Meza
    Carlos Manzano
*/
#include <cstdlib>
#include <iostream>
#include <conio.h>
#include <stdlib.h>
#include <iomanip>
#include <string>
#include <stdio.h>
#include <math.h>

#define f1(x)(1.0 * pow(x,4) - 8.0 * pow(x,3) - 35.0 * pow(x , 2)  + 450.0 * x - 1001)  //secante
#define fb(x)( ( 3.0 * 3.1416 * pow(x,2) ) - ( (3.1416 * pow(x,2)) / 3) - 30 ) //Biseccion  //ERROR
#define f2(x)(10.0 - 20.0*(exp(-0.15 * x) - exp(-0.5 * x))-5)
#define fp(x)(3.0 * pow(0.6065,x)*(pow(1.419,x) - 3.33))
#define f3(x)(1.0 * pow(x,4) - 5.0 * pow(x,3) + 10.0 * pow(x,2) - 10.0 * x + 4.0)

using namespace std;

void leemat(int , float [30][30]);
void escmat(int, float [30]);
float gaussjor( int , float [30][30], float [30] );
void PideDatos(int *Dat, int *Ord, float Val[][102]);
float Potencia(int n, float Num);
void PreparaSistema(int Ord, int Dat, float Sist[][102], float Val[][102]);
void ResuelveGauss(int Dim, float Sist[][102]);
void EscribeDatos(int Dim, float Sist[][102]);
void leevec(int, float [100] , char[100]);
void DifDivNewton(int, float [100], float [100], float [100]);
void escribirMatriz(int , int , float [50][50], char[10]);
void gaussjordan(int , float [50][50], float [50][50], float [50]);
void escvec(int, float [50], char[10]);
void escpol(int, float [50], char [10]);
void intlag(int, float[50], float[50][50]);
double fbi(double x);
double biseccion ( double a, double b, double tol, int maxlter);

int Menu(){
    int iOpcion;

    cout << endl;
    cout << "Seleccione el numero del metodo que desee: " << endl;
    cout << "   1.- Trapezoidal" << endl;
    cout << "   2.- Minimos cuadrados" << endl;
    cout << "   3.- Simpson 1/3" << endl;
    cout << "   4- Simpson 3/8" << endl;
    cout << "   5.- Secante" << endl;
    cout << "   6.- Newton-Raphson" << endl;
    cout << "   7.- Biseccion" << endl;
    cout << "   8.- Lagrange" << endl;
    cout << "   9.- Eliminacion Gausseana" << endl;
    cout << "   10.- Newton (diferencias divididas) " << endl;
    cout << "   11.- Montante" << endl;
    cout << "   12.- Salir" << endl;


    cout << endl << "Opcion: ";
    cin >> iOpcion;
    cout << endl;

    return iOpcion;
}

void monty(){
	int n, x, y,w;
	cout<<"Ingrese el grado de la matriz: ";
	cin>> n;

	double A[n][n];
	double B[n][n];
	double mata[n];
	double matb[n];

	for (x = 0; x < n; x++){
			cout<< "Ingrese los elementos del renglon " << x+1 <<endl;
		for(y = 0; y < n; y++){
			cin>> A[x][y];
			B[x][y]=0;
		}
		matb[x]=0;
	}
	cout<<"Ingrese la columna de resultados:" << endl;
	for (x = 0; x < n; x++){
		cin >> mata[x];
		matb[x] = mata[x];
	}
	cout<<"Matriz seleccionada" <<endl;
	for (x = 0; x < n; x++){ //MUESTRA LA MATRIZ INGRESADA
		for(y = 0; y < n; y++){
			cout<<A[x][y]<<"\t";
		}
		cout<<endl;
	}
	for(x = 0; x < n ; x++){
		for (y = 0; y < n ; y++){
			B[x][y]=A[x][y];
		}
		for (y = 0; y < n ; y++){
			if (y == x){
			}
			else {
				B[y][x]=0;
			}
		}
		for (y = 0; y < n; y++){
			B[y][y]=A[x][x];
		}
		for (y = 0; y < n; y++){
			if (y == x){
			}
			else {
				matb[y] = (mata[y]* A[x][x]) - (A[y][x] * mata[x]);
				if (x == 0){
				}
				else {
					matb[y]= matb[y] / A[0][0];
				}
			}
			for (w = x+1 ; w<n ; w++){
			if (y == x){
			}
			else {
				B[y][w]= (A[y][w]* A[x][x]) - (A[y][x] * A[x][w]);
				if (x == 0){
				}
				else {
					B[y][w]=B[y][w] / A[0][0];
				}
			}
		}
		}
		for (y = 0; y < n ; y++){
			for (w = 0; w<n; w++){
				A[y][w]=B[y][w];
			}
			mata[y]=matb[y];
		}
	}

cout<<"Los valores de x son : ";
    for(x = 0; x<n ; x++){
        cout<<"X"<<x<<" = "<<matb[x] / B[x][x] <<endl;
    }
}

float f(float(x))
{
    return (pow(x,3)+pow(x,2)-(4*x)-5);
}


float g(float(x)){
    return (3*pow(x,2)+2*x-4);
}


float h(float(x)){
    return (6*x+4);
}

double funcion(double dx){
       double dFuncion;
       dFuncion= 0.2 + 25.0 * dx - 200.0 * pow(dx ,2) + 675.0 * pow(dx , 3) - 900.0 * pow(dx , 4) + 400.0 * pow(dx, 5);
       return dFuncion;
}


double simpson ( double dArea, double dS1, double dS2, double dx, double dN, double dA, double dB, double di, double dH){
    cout << "Ingrese el numero de Trapecios: ";
    cin >> dN;
    cout << "Ingrese el limite inferior: ";
    cin >> dA;
    cout << "Ingrese el limite superior ";
    cin >> dB;

    dS1 = 0;
    dS2 = 0;
    dx = dA;
    dH = (dB - dA) / dN;

    cout<< "h: " << dH << endl;

    if(dN == 2) {
        dx = dx + dH;
        dS1 = dS1 + funcion(dx);
        dArea = (dH / 3) * (funcion(dA) + funcion(dB) + 4 * dS1 + 2 * dS2);
        cout << "El area es : " << dArea;
        cout << endl;
    }
    else
        {
            for(di = 1; di <= ((dN / 2 ) - 1); di++) {
                dx = dx + dH;
                dS1 = dS1 + funcion(dx);
                dx = dx + dH;
                dS2 = dS2 + funcion(dx);
            }
    dx = dx + dH;
    dS1 = dS1 + funcion(dx);
    dArea =((funcion(dA) + funcion(dB) + 4 * dS1 + 2 * dS2) *(dH / 3));
    cout << "El area es : " << dArea;

      }

}

//SECANTE
double secante(double dX0, double dX1, double des, int iIter){
    double dX2, dea, dY0, dY1, dPend;

    cout << "iteracion" << setw(10) << "xi-1" << setw(10) <<
    "xi" << setw(10) << "xi+1" << setw(10) << "Error" << setw(10)<<
    "f(xi-1)"<< setw(10) << "f(xi)" << setw(10) << "pend" << endl;

    for(int iCont = 1; iCont <= iIter; iCont++){
        dY0 = f1(dX0);
        dY1 = f1(dX1);
        dPend = (dY1 - dY0)/(dX1 - dX0);
        dX2 = dX1 - (dY1/dPend);
        dea = fabs((dX2 - dX1)/(dX2));

        cout << setw(5) << iCont << setw(15)<< dX0 << setw(10) << dX1 << setw(10)<< dX2 << setw(10)<<
        dea<< setw(10) << dY0 << setw(10)<< dY1<< setw(10) << dPend << endl;

        dX0 = dX1;
        dX1 = dX2;

        if(dea < des) {
            cout << "El metodo converge a las " << iCont << " iteraciones" << endl;
            return dX2;
            break;
        }

    }
    if( dea > des)
        cout << "El metodo no converge para las iteraciones especificadas" << endl;
        return 0;

}

//BISECCION
 double fbi(double x){
        return  (( 3.0 * 3.1416 * pow(x,2) ) - ( (3.1416 * pow(x,2)) / 3) - 30 );
 }
 double biseccion(double a, double b, double tol, int maxlter){
        double c;
        int nolter = 0;
        do {
            c = (a+b)/2;
            if(fbi(a)*fbi(c)<0)
            {
               b = c;
            }
            else
               a = c;

            cout << nolter << "\t" << a << "\t" << b << "\t" << c << "\t" <<f(c)<<endl;
            nolter++;
         }
         while((fabs(fbi(c))>tol)&&(nolter<maxlter));
         return c;
 }

//NEWTON - RAPHSON
double NewtonRaphson(double LimiteInferior, double errorMenor, int numeroIteracciones ){
    double x2, ea, y1, yp1;
    cout << "iteracion " << setw(10) << "xi" << setw(10) << "xi+1" << setw(10) << "Ea" << setw(10) << "f(x1)" << setw(10) << "f'(x1)" << endl;
    for(int cont = 1; cont <= numeroIteracciones; cont++)
    {
        y1 = f2(LimiteInferior);
        yp1 = fp(LimiteInferior);
        if(yp1 == 0)
        {
            cout << "f'(xi) = 0" << endl;
            break;
        }
        x2 = LimiteInferior - ( y1 / yp1 );
        ea = fabs(( x2 - LimiteInferior)/x2);

        cout << setw(5) << cont << setw(15) << LimiteInferior << setw(10) << x2 << setw(10) << ea << setw(10) << y1 << setw(10) << yp1 << endl;

        LimiteInferior = x2;

        if(ea < errorMenor)
        {
            cout << "El metodo coincide a " << cont << " iteraciones." << endl;
            return x2;
            break;
        }
    }
    if(ea > errorMenor) {
        cout << "El metodo no coincide para las " << numeroIteracciones << " iteraciones indicadas" << endl;
    }
    return 0;
}

float gaussjor(int n, float a[30][30], float x[30]){
    int iA, j, k;
    float aux;
    for( iA = 1; iA <= n; iA++ ){
        aux = a[iA][iA];
        for ( j = iA; j <= n + 1; j++ ){
            a[iA][j] = a[iA][j] / aux;
        }
        for ( j = 1; j <= n; j++ ){
            if( j != iA ){
                aux = a[j][iA];
                for (k = iA; k <= n + 1; k++){
                    a[j][k] = ( a[j][k] ) - ( ( a[iA][k] ) * aux );
                }
            }
        }
    }

    for ( iA = 1; iA <= n; iA++ )
        x[iA] = a[iA][n + 1];
}

void leemat(int n, float a[30][30] ){
    int i, j;
    for (i = 1; i <= n; i++ ){
        for( j = 1; j <= n + 1; j++ ){
            cout << "a[" << i << "][" << j << "] = ";
            cin >> a[i][j];
        }
    }
}

void escmat( int n, float x[30] ){
    int i;
    for (i = 1; i <= n; i++)
        cout << "x[" << i << "] = " << x[i] << endl;
}

double Trapezoidal(){
    int N;
    double h,a,b,i,suma,F,xi,t;

    cout<<"METODO DE INTEGRACION TRAPEZOIDAL"<<endl<<endl<<endl;
    cout<<"LA FUNCION A INTEGRAR ES: SEN(X)/(1+X^4)^0.5 "<<endl<<endl<<endl;
    cout<<"INGRESE EL LIMITE INFERIOR DE LA INTEGRAL--(A): ";
    cin>>a;//LIMITE INFERIODE LA INTEGRACION
    cout<<"INGRESE EL LIMITE SUPERIOR DE LA INTEGRAL--(B): ";
    cin>>b;//LIMITESUPERIO DE LA INTEGRACION
    cout<<"INGRESE N: ";
    cin>>N;//ES EL NUMERO DE LAS ITERACIONES

    h=(b-a)/N;
    i=0;
    suma=0;

    for(i=0;i<=N;i++){
        xi=a+i*h;
        F=sin(xi)/ pow( (1+ xi*xi*xi*xi) , 0.5);
        suma=suma+F;
        t=suma*h;
    }

    cout<<"los valores son : "<<endl ;
    cout<<"a= "<<a<<endl; // endl hace salto de linea
    cout<<"b= "<<b<<endl;
    cout<<"h= "<<h<<endl<<endl;
    cout<<"suma es: "<<suma<<endl<<endl<<endl;
    cout<<" LA INTEGRAL ES : "<<t<<endl;

}

void PideDatos(int *Dat, int *Ord,float Val[][102])
{
    int A,B;

    printf("\n\n\n METODO DE MINIMOS CUADRADOS.\n\n");
    printf("\n Introduce el numero de datos (Puntos): ");scanf("%d",&*Dat);
    printf("\n\n\n Introce los valores de cada punto\n");

    for(A=1;A<=*Dat;A++)
    {
        printf(" -Valores del Punto %d:\n",A);
        printf("   X%d: ",A); scanf("%f",&Val[0][A]);
        printf("   Y%d: ",A); scanf("%f",&Val[1][A]);
    }
    printf("\n\n\n Introduce el orden del polinomio: "); scanf("%d",&B);
    *Ord=B+1;
}


float Potencia(int n, float Num)
{
    int A;
    float res;

    res=1;
    for(A=1;A<=n;A++) res=res*Num;
    return res;
}



void PreparaSistema(int Ord, int Dat, float Sist[][102], float Val[][102])
{
    int A,B,C,Exp;
    float suma,termino;

    for(A=1;A<=Ord;A++)    for(B=1;B<=Ord;B++)
    {
        suma=0;
        Exp=A+B-2;

        for(C=1;C<=Dat;C++)
        {
            termino=Val[0][C];
            suma=suma+Potencia(Exp,termino);
        }
        Sist[A][B]=suma;
    }
    for(A=1;A<=Ord;A++)
    {
        suma=0;
        Exp=A-1;

        for(C=1;C<=Dat;C++)
        {
            termino=Val[0][C];
            suma=suma+Val[1][C]*Potencia(Exp,termino);
        }
        Sist[A][Ord+1]=suma;
    }
}


void ResuelveGauss(int Dim, float Sist[][102])
{
    int NoCero,Col,C1,C2,A;
    float Pivote,V1;

    for(Col=1;Col<=Dim;Col++){
        NoCero=0;A=Col;
        while(NoCero==0){
            if(Sist[A][Col]!=0){
                NoCero=1;}
            else A++;}
        Pivote=Sist[A][Col];
        for(C1=1;C1<=(Dim+1);C1++){
            V1=Sist[A][C1];
            Sist[A][C1]=Sist[Col][C1];
            Sist[Col][C1]=V1/Pivote;}
        for(C2=Col+1;C2<=Dim;C2++){
            V1=Sist[C2][Col];
            for(C1=Col;C1<=(Dim+1);C1++){
                Sist[C2][C1]=Sist[C2][C1]-V1*Sist[Col][C1];}
        }}

    for(Col=Dim;Col>=1;Col--) for(C1=(Col-1);C1>=1;C1--){
        Sist[C1][Dim+1]=Sist[C1][Dim+1]-Sist[C1][Col]*Sist[Col][Dim+1];
        Sist[C1][Col]=0;
    }
}

void EscribeDatos(int Dim, float Sist[][102])
{
    int A,B;
    printf("\n\n");
    for(A=1;A<=Dim;A++){
        for(B=1;B<=(Dim+1);B++){
            printf("%7.2f",Sist[A][B]);
            if(B==Dim) printf("   |");}
        printf("\n");
    }
    printf("\n\n");
}

void leevec(int iN, float iX[100], char cNom[10])
{
    for(int iC = 1; iC <= iN ; iC++)
    {
        cout << cNom << "[" << iC << "] : ";
        cin >> iX[iC];
    }
}

void escribirMatriz(int iN , int iM, float dA[50][50], char cNom[10])
{
    cout << "Matriz" << cNom << endl;
    for (int iC = 1; iC <= iN; iC++)
    {
        cout << endl;
        for(int iR = 1; iR <= iM; iR++)
        {
            cout << setw(10) << dA[iC][iR];


        }
    }
    cout << endl;
}

void escvec(int iN, float dX[50], char cNom[10])
{
    for(int iC = 1; iC <= iN ; iC++)
    {
        cout << endl;
        cout << cNom << "[" << iC << "] : " << dX[iC];
    }
    cout << endl;
}

void DifDivNewton(int iN, float dX[100], float dY[100], float dA[100])
{
    float DD[50][50], dP[50][50], dVX[50][50], dB[50], P[50][50];

    for(int iC = 1; iC <= iN ;iC++)
    {
        DD[iC][1] = iC - 1;
        DD[iC][2] = dX[iC];
        DD[iC][3] = dY[iC];
    }
    for(int iR = 1; iR <= iN; iR++ )
    {
        for(int iC = 1; iC <= iN - iR; iC++)
        {
            DD[iC][iR+3] = ((DD[iC+1][2+iR])- (DD[iC][2 + iR]))/((DD[iC + iR][2]) -(DD[iC][2]));
        }
    }
    escribirMatriz(iN, iN+2,DD, "DD");
    for (int iC = 1; iC <= iN; iC++)
    {
        dB[iC] = DD[1][iC+2];
    }

    cout << "Valores de B" << endl;
    escvec(iN, dB, "dB");
    for(int iC = 1; iC <= iN ; iC++)
    {
        for(int iR = 1; iR <= iN; iR++)
        {
            P[iC][iR] = 0;
        }

    }
    for(int iC = 1; iC <= iN ; iC++)
    {
        P[iC][1] = 1;
    }
    for(int iC = 2; iC <= iN; iC++)
    {
        for(int iR = 2; iR <= iC; iR++)
        {
            P[iC][iR] = ((P[iC -1][iR - 1])* (-1 * dX[iC -1 ]))+ P[iC - 1][iR];
        }
    }

    escribirMatriz(iN, iN, P, "Polinomios");


    for(int iC = 1; iC <= iN; iC++)
    {
        for(int iR = 1; iR <= iN; iR++)
        {
            P[iC][iR] = P[iC][iR] * dB[iC];
        }
    }
    escribirMatriz(iN, iN , P , "Polinomios");

    for(int iR = 1; iR <= iN; iR++)
    {
        dA[iR -1] = 0;
        for(int iC = iR; iC <= iN ; iC++)
        {
            dA[iR -1] = dA[iR - 1] + P[iC][iC+1 - iR];
        }
    }
    escpol(iN - 1, dA, "P(x)" );

}

void escpol(int iN, float dA[50], char cNom[10])
{
    cout << "Polinomio" << cNom << " = ";
    for(int iC = 0; iC <= iN; iC++)
    {
        cout << dA[iC] << "X^" << iC;
            if (iC < iN)
            {
                cout << "+ ";
            }
    }
    cout << endl;
    }


void intlag(int n, float a[50], float fx[50][50] ){
    int i, j, k, cont, h, l, t;
    float P[50][50], constP[50][50], x[50], aux[50], PL[50];

    cont = 1;
    for ( i = 0; i <= 50; i++ ){
        aux[i] = 1;
    }

    for ( j = 0; j < n; j++ ){
        for (k = 0; k < n; k++ ){
            if ( j != k )
                aux[j] = aux[j] * ( a[j]- a[k] );
            fx[j][2] = fx[j][1] / aux[j];
        }
    }

    for ( i = 0; i < n; i++ )
        a[i] = a[i] * -1;

    for ( l = 0; l <= n - 1; l++ ){
        t = 0;
        for ( j = 0; j < n; j++ ){
            if (j != l){
                x[t] = a[j];
                t++;
            }
        }

        P[0][0] = 1;
        for (i = 0;  i < n; i++ ){

            for ( j = 1; j < n; j++ ){
                if (i == j)
                    P[i][j] = P[i - 1][j - 1] * x[i - 1];
            }

            for ( j = 1; j <= n; j++ ){
                if (i == 0)
                    P[j][i] = P[0][0];
            }

            for ( j = 1; j < n - 1; j++ ){
                if (i == 1)
                    P[j + 1][i] = P[j][i] + x[j];
            }

            for ( j = 2; j < n; j++ ){
                for( h = j + 1; h < n; h++ ){
                    if (i == j)
                        P[h][j] = P[h - 1][j] + (P[h - 1][j - 1] * x[h - 1]);
                }
            }
        }

        for (j = 0; j < n; j++ )
            constP[0][j] = fx[l][2] * P[n - 1][j];

        for (i = 0; i < n; i++ ){
            for (j = 0; j < n; j++){
                PL[i] = PL[i] + constP[j][i];
            }
        }
        cont++;
    }

    cout << endl << endl;;
    cout << "El polinomio de Lagrange: " << endl << endl;
    cout << "(" <<PL[0]-1 << "x^" << n-1-i << ") + ";
    for (i = 1; i < n; i++ )
        cout << "(" <<PL[i] << "x^" << n-1-i << ") + ";
}

int main()
{
    int iOpcion, numeroIteracciones;
    double errorMenor, LimiteInferior, LimiteSuperior, Raiz, dX0, dX1;
    int n;
    float a[30][30], x[30];

    do {
        iOpcion = Menu();
        if (iOpcion == 1){
            Trapezoidal();
        }
        else if (iOpcion == 2){
            int Datos,Orden,C;

            float Valores[2][102],Sistema[102][102];
            PideDatos(&Datos,&Orden,Valores);
            PreparaSistema(Orden,Datos,Sistema,Valores);
            printf("\n\n El sistema a resolver es el siguiente:");
            EscribeDatos(Orden,Sistema);
            ResuelveGauss(Orden,Sistema);
            printf("\n\n El sistema resuelto:");
            EscribeDatos(Orden,Sistema);
            printf("\n\n La Ecuacion del Polinomio ajustado por minimos Cuadrados\n\n:");
            for(C=1;C<=Orden;C++)
                printf(" + (%f)X^%d",Sistema[C][Orden+1],C-1);
        }
        else if (iOpcion == 3){
            double  dArea, dS1, dS2, dx , dN, dA, dB, di , dH;
            simpson ( dArea,dS1,dS2,dx,dN,dA,dB,di,dH);
        }
        else if (iOpcion == 4){
            int i;
            float a, b, d, n, I = 0, J = 0, A, K = 0, E = 0;
            cout << " Formula f(x)= x^3 + 2x^2 - 4x - 5 " << endl;
            cout << "Limite inferior: ";
            cin >> a;
            cout << endl;
            cout << "Limite superior ";
            cin >> b;
            cout << endl;
            cout << "Numero de intervalos : ";
            cin >> n;
            cout << endl;
            d = (b - a) / n;

            for(i = 1; i < n; i++) {
                I = I + f(a + (i * d));
            }
            for(i = 3; i < n - 1; i++){
                if((i % 3) == 0){
                    J = J + f(a + (i * d));
                }
            }

            A = 3 *(d / 8) * (f(a) + (3 * I) - J + f(b));

            cout << "El valor de la integral es: ";
            cout << A <<endl;

            E = -(d * d * d * d * d * 6 * 3 / 80);
//
            cout << "El error total es : ";
            cout << E << endl;
        }
        else
            if (iOpcion == 5){
            int iIter;
            double des, dX0, dX1, dRaiz;

            cout << "Punto inicial (xi-1)" << endl;
            cin >> dX0;
            cout << "Punto inicial 2 (xi)" << endl;
            cin >> dX1;
            cout << "Error"<< endl;
            cin >> des;
            cout << "Numero de iteraciones" << endl;
            cin >> iIter;
            cout.precision(4);

            dRaiz = secante(dX0, dX1, des, iIter);
            cout << "La raiz es: " << setprecision(10) << dRaiz << endl;
            cout << endl;
        }
        else if (iOpcion == 6){
            cout << endl;
            cout << "Punto Inicial (Xi)= ";
            cin >> LimiteInferior;
            cout << endl;
            cout << "Error = ";
            cin >> errorMenor;
            cout << endl;
            cout << "Numero de iteracciones (especifique si el metodo termina antes) = ";
            cin >> numeroIteracciones;
            cout << endl;
            cout.precision(4);
            Raiz = NewtonRaphson( LimiteInferior, errorMenor,  numeroIteracciones );
            cout << "La raiz es: " <<setprecision(10) << Raiz;
            cout << endl;
        }
        else if (iOpcion == 7 ){
            double a, b, tol, raiz;
            int maxlter;

            cout << "Metodo para Biseccion" << endl << endl;
            cout << "Limite Inferior:  ";
            cin >> a;
            cout << endl;
            cout << "Limite superior:  ";
            cin >> b;
            cout << endl;
            cout << "Error:  ";
            cin >> tol;
            cout << endl;
            cout << "Iteraccion:  ";
            cin >> maxlter;
            cout << endl;
            raiz = biseccion(a,b,tol,maxlter);
            cout << endl;
            cout << "La raiz es: " << raiz << "%" << endl;
        }
        else if ( iOpcion == 8 ){
            int n, i;
            float x[50], fx[50][50];

            cout << endl;
            cout << "Metodo de interpolacion de lagrange" << endl;
            cout << "Numero de datos de la lista" ;
            cin >> n;
            cout << "Introduce los datos de x[i]: " << endl;
            for ( i = 0; i <= n - 1; i++){
                cout << endl;
                cout << "x[" << i <<"]= ";
                cin >> x[i];
            }
            cout << endl;
            cout << "Introduce los valores de f(x)" << endl;
            for( i = 0; i <= n - 1; i++ ){
                cout << endl;
                cout << "fx(" <<i <<")= ";
                cin >> fx[i][1];
            }
            intlag(n, x ,fx);
        }
        else if ( iOpcion == 9 ){
            cout << "Numero de ecuaciones: ";
            cin >> n;
            cout << "Introduce los valores de la matriz: ";
            leemat(n, a);
            gaussjor(n, a, x);
            cout << "Soluccion " << endl;
            escmat(n, x);
        }
        else if ( iOpcion == 10 ){
            int iM, iN;
            float dX[100], dY[100], dA[100];

            cout << "Cuantos puntos desea : " << endl;
            cin >> iN;
            cout << "Vector de x" << endl;
            leevec(iN, dX, "X");
            cout << "Vector de y" << endl;
            leevec(iN, dY, "Y");

            DifDivNewton(iN, dX, dY, dA);

        }
        else if ( iOpcion == 11 ){
            monty();
        }
        else if ( iOpcion == 12 ){
            return 0;
        }
    }
    while ( iOpcion < 12 || iOpcion > 0 );

    return 0;
}
