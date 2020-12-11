// 714HW1.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include "pch.h"
#include <iostream>
#include <stdio.h>
#include <math.h>       
#include "Grid.h"


#define PI 3.14159265

int main()
{
	//ErrorOfLinearInterpolation();
	

	/*ErrorVsGridSpacing2DWave();

	ErrorVsGridSpacing2DWaveSpectral();*/
	
	ErrorVsGridSpacing2DCooling();
	ErrorVsTimeSteps2DCooling();

	//int mRef = 31; //255

	//Grid oReference(mRef);
	//oReference.ExactSolution_Cooling();
	//oReference.OutputResultAsVTK();
	
	//oReference.Solve2DCoolingByADI();
	//oReference.Solve2DCoolingByLOD();
	//oReference.Solve2DCoolingByCrankNicolson();

	//ErrorVsGridSpacing1DInvisBurgers();

	//oReference.Solve1DInvisBurgers_nonConServative();
	//oReference.Solve1DInvisBurgers_ConServative();
	//oReference.Solve1DInvisBurgers_Godunov();

	

	//ErrorVsGridSpacing();
	
	//// multi grid
	//int mRef = 63; //255
	//Grid oReference(mRef);
	//oReference.SolveByMultiGrid(1e5);
	//oReference.OutputResultAsVTK();

	getchar();
}

void ErrorOfLinearInterpolation()
{
	int iCheckN = 50;

	double dH;
	double dH_check ;

	double df_0, df_1;
	double dx_0, dx_1;
	double dx;
	double dInterpolation;

	double dError;
	double dMaxError ;



	for (int iN = 1; iN <= 200; iN++)
	{
		//int iN=2;
		

		dH = 1. / iN;
		dH_check = dH / iCheckN;
		dMaxError = 1.e-16;


		for (int i = 1; i <= iN; i++)
		{
			dx_0 = (i - 1) * dH;
			dx_1 = i * dH;

			df_0 = FunctionX(dx_0);
			df_1 = FunctionX(dx_1);

			for (int j = 1; j <= iCheckN; j++)
			{
				dx = dx_0 + j * dH_check;
				dInterpolation = df_0 + (dx - dx_0) / (dx_1 - dx_0) * (df_1 - df_0);

				dError = fabs(FunctionX(dx) - dInterpolation);

				//printf("f=%lg , dInterpolation =%lg \n", FunctionX(dx), dInterpolation);

				if (dError > dMaxError) dMaxError = dError;
			}
		}

		printf("N=%d , Error =%lg \n", iN, dMaxError);

		if (dMaxError < 0.01) break;
	}
}

const int iB = 8;

double FunctionX(double x)
{
	//return exp(-400. * (x - 0.5)* (x - 0.5));

	return sin(iB*PI*x);
}

double FunctionXSpectral(double x)
{
	
	//return exp(-100. * (x)* (x));

	return sin(iB*PI*(x/2. + 0.5));
}


void ErrorVsGridSpacing2DWaveSpectral()
{
	int mRef = 63; //255
	int iRatio = 2;

	int mGrid_1 = (mRef + 1) / iRatio - 1;

	Grid oReference(mRef);
	oReference.SolveTwoDWavebySpectral();

	for (int iL = 0; iL <= 3; iL++)
	{
		Grid oGrid_1(mGrid_1);
		oGrid_1.SolveTwoDWavebySpectral();


		// compara error
		double MaxError = 1.e-30;
		double Error;

		for (int i = 0; i <= mGrid_1 + 1; i++)
			for (int j = 0; j <= mGrid_1 + 1; j++)
			{
				{
					Error = fabs(oGrid_1.m_pUij[i][j] - oReference.m_pUij[i*iRatio][j*iRatio]);
					if (Error > MaxError) MaxError = Error;
				}
			}

		printf("Spectral, Error beteen grid =%d and %d is %lg \n", mGrid_1, mRef, MaxError);


		iRatio *= 2;
		mGrid_1 = (mRef + 1) / iRatio - 1;
	}
}

void ErrorVsGridSpacing2DWave()
{
	int mRef = 511; //255 //1023
	int iRatio = 2;

	int mGrid_1 = (mRef + 1) / iRatio - 1;

	Grid oReference(mRef);
	oReference.SolveTwoDWave();

	for (int iL = 0; iL <= 3; iL++)
	{
		Grid oGrid_1(mGrid_1);
		oGrid_1.SolveTwoDWave();
		

		// compara error
		double MaxError = 1.e-30;
		double Error;

		for (int i = 0; i <= mGrid_1 + 1; i++)
			for (int j = 0; j <= mGrid_1 + 1; j++)
			{
				{
					Error = fabs(oGrid_1.m_pUij[i][j] - oReference.m_pUij[i*iRatio][j*iRatio]);
					if (Error > MaxError) MaxError = Error;
				}
			}

		printf("FDM, Error beteen grid =%d and %d is %lg \n", mGrid_1, mRef, MaxError);

		iRatio *= 2;
		mGrid_1 = (mRef + 1) / iRatio - 1;
	}

}

void ErrorVsTimeSteps2DCooling()
{
	double dt = 1.e-2; // 1.e-4

	int mRef = 127; //255
	int iRatio = 2;
	int iLoopMax = 1e6;

	int mGrid_1 = mRef;



	for (int iL = 0; iL <= 3; iL++)
	{

		

		Grid oGrid_1(mGrid_1);
		//oGrid_1.Solve2DCoolingByADI(dt);

		//oGrid_1.Solve2DCoolingByLOD(dt);
		oGrid_1.Solve2DCoolingByCrankNicolson(dt);
		//oGrid_1.Solve2DCoolingByImplicitEuler(dt);

		Grid oAnalytical(mGrid_1);
		oAnalytical.ExactSolution_Cooling();

		// compara error
		double MaxError = 1.e-30;
		double Error;

		for (int i = 0; i <= mGrid_1 + 1; i++)
			for (int j = 0; j <= mGrid_1 + 1; j++)
			{
				{
					Error = fabs(oGrid_1.m_pUij[i][j] - oAnalytical.m_pUij[i][j]);
					if (Error > MaxError) MaxError = Error;
				}
			}

		printf("Grid =%d, Error for dt =%lg is %lg , nTDM=%d\n", mGrid_1,dt, MaxError, oGrid_1.m_inTDM);

		
		dt = dt / iRatio;
	}

}

void ErrorVsGridSpacing2DCooling()
{
	double dt = 1.e-6; // 1.e-4

	int mRef = 255; //255
	int iRatio = 2;
	int iLoopMax = 1e6;

	int mGrid_1 = (mRef + 1) / iRatio - 1;

	

	for (int iL = 0; iL <= 2; iL++)
	{
		
		dt = 0.01*1. / (mGrid_1 + 1); //0.01 h
		//dt = 1. / (mGrid_1 + 1) / (mGrid_1 + 1);  //h^2
		//dt = 0.01 / (mGrid_1 + 1) / (mGrid_1 + 1);

		Grid oGrid_1(mGrid_1);
		oGrid_1.Solve2DCoolingByADI(dt);

		//oGrid_1.Solve2DCoolingByLOD(dt);
		//oGrid_1.Solve2DCoolingByCrankNicolson(dt);
		//oGrid_1.Solve2DCoolingByImplicitEuler(dt);
		
       Grid oAnalytical(mGrid_1);
	   oAnalytical.ExactSolution_Cooling();

		// compara error
		double MaxError = 1.e-30;
		double Error;

		for (int i = 0; i <= mGrid_1 + 1; i++)
			for (int j = 0; j <= mGrid_1 + 1; j++)
			{
				{
					Error = fabs(oGrid_1.m_pUij[i][j] - oAnalytical.m_pUij[i][j]);
					if (Error > MaxError) MaxError = Error;
				}
			}

		printf("Error for grid =%d is %lg , nTDM=%d\n", mGrid_1, MaxError , oGrid_1.m_inTDM);

		iRatio *= 2;
		mGrid_1 = (mRef + 1) / iRatio - 1;
	}

}


void ErrorVsGridSpacing1DInvisBurgers()
{
	int mRef = 255; //255
	int iRatio = 2;
	int iLoopMax = 1e6;

	int mGrid_1 = (mRef + 1) / iRatio - 1;

	Grid oReference(mRef);
	oReference.Solve1DInvisBurgers_nonConServative();
	
	for (int iL = 0; iL <= 2; iL++)
	{
		Grid oGrid_1(mGrid_1);
		oGrid_1.Solve1DInvisBurgers_nonConServative();
		

		// compara error
		double MaxError = 1.e-30;
		double Error;

		
			for (int j = 0; j <= mGrid_1 + 1; j++)
			{
				
					Error = fabs(oGrid_1.m_pA[j] - oReference.m_pA[j*iRatio]);
					if (Error > MaxError) MaxError = Error;
				
			}

		printf("Error beteen grid =%d and %d is %lg \n", mGrid_1, mRef, MaxError);

		iRatio *= 2;
		mGrid_1 = (mRef + 1) / iRatio - 1;
	}

}

void ErrorVsGridSpacing()
{
	int mRef = 255; //255
	int iRatio = 4;
	int iLoopMax = 1e6;

	int mGrid_1 = (mRef + 1) / iRatio - 1;

	Grid oReference(mRef);
	oReference.SolveByJacobi(iLoopMax, 1.);
	oReference.OutputResultAsVTK();

	for (int iL = 0; iL <= 2; iL++)
	{
		Grid oGrid_1(mGrid_1);
		oGrid_1.SolveByJacobi(iLoopMax, 1.);
		oGrid_1.OutputResultAsVTK();


		// compara error
		double MaxError = 1.e-30;
		double Error;

		for (int i = 0; i <= mGrid_1 + 1; i++)
			for (int j = 0; j <= mGrid_1 + 1; j++)
			{
				{
					Error = fabs(oGrid_1.m_pUij[i][j] - oReference.m_pUij[i*iRatio][j*iRatio]);
					if (Error > MaxError) MaxError = Error;
				}
			}

		printf("Error beteen grid =%d and %d is %lg \n", mGrid_1, mRef, MaxError);

		iRatio *= 2;
		mGrid_1 = (mRef + 1) / iRatio - 1;
	}

}

// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu

// Tips for Getting Started: 
//   1. Use the Solution Explorer window to add/manage files
//   2. Use the Team Explorer window to connect to source control
//   3. Use the Output window to see build output and other messages
//   4. Use the Error List window to view errors
//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file
