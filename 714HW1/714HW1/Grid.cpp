#include "pch.h"
#include "Grid.h"
#include <iostream>

#define PI 3.14159265

Grid::Grid()
{
}

Grid::Grid(int m)
{
	m_iMeshM = m;

	m_Alpha = 1.;
	m_TargetTime = 0.01;
	
	m_inTDM = 0;
	m_isExact = 0;

	m_dH = 1. / (m_iMeshM + 1);
	m_iNtotal = (m_iMeshM + 2)*(m_iMeshM + 2);

	m_iChebyN = m_iMeshM + 1;

	m_iT = 0;

	m_pUij = new double *[m_iMeshM + 2];   //m+2*m+2    0~m+1
	for (int i = 0; i <= m_iMeshM + 1; i++) m_pUij[i] = new double[m_iMeshM + 2];

	m_pU0ij = new double *[m_iMeshM + 2];
	for (int i = 0; i <= m_iMeshM + 1; i++) m_pU0ij[i] = new double[m_iMeshM + 2];

	m_pRij = new double *[m_iMeshM + 2];
	for (int i = 0; i <= m_iMeshM + 1; i++) m_pRij[i] = new double[m_iMeshM + 2];


	m_pU00ij = new double *[m_iMeshM + 2];
	for (int i = 0; i <= m_iMeshM + 1; i++) m_pU00ij[i] = new double[m_iMeshM + 2];



	for (int i = 0; i <= m_iMeshM + 1; i++)
		for (int j = 0; j <= m_iMeshM + 1; j++)
		{
			{
				m_pUij[i][j] = 0.;
				m_pU0ij[i][j] = 0.;
				m_pRij[i][j] = 0.;

				m_pU00ij[i][j] = 0.;
			}
		}


		m_pA = new double[m_iMeshM + 2];  // 0~m+1
	    m_pB = new double[m_iMeshM + 2];
		m_pC = new double[m_iMeshM + 2];
		m_pD = new double[m_iMeshM + 2];

		for (int i = 0; i <= m_iMeshM + 1; i++)
		
		{
			m_pA[i] = 0.;
			m_pB[i] = 0.;
			m_pC[i] = 0.;
			m_pD[i] = 0.;
			
		}
		


}

Grid::~Grid()
{
	if (m_pUij != NULL)
	{
		for (int i = 0; i <= m_iMeshM + 1; i++)
		{
			delete[] m_pUij[i];
			m_pUij[i] = NULL;
		}
		delete[] m_pUij;
		m_pUij = NULL;
	}

	if (m_pU0ij != NULL)
	{
		for (int i = 0; i <= m_iMeshM + 1; i++)
		{
			delete[] m_pU0ij[i];
			m_pU0ij[i] = NULL;
		}
		delete[] m_pU0ij;
		m_pU0ij = NULL;
	}

	if (m_pRij != NULL)
	{
		for (int i = 0; i <= m_iMeshM + 1; i++)
		{
			delete[] m_pRij[i];
			m_pRij[i] = NULL;
		}
		delete[] m_pRij;
		m_pRij = NULL;
	}

	if (m_pU00ij != NULL)
	{
		for (int i = 0; i <= m_iMeshM + 1; i++)
		{
			delete[] m_pU00ij[i];
			m_pU00ij[i] = NULL;
		}
		delete[] m_pU00ij;
		m_pU00ij = NULL;
	}
}

void Grid::TwoDWaveStep2()
{

	int iCount;
	double Error, MaxError;
	double x, y;
	double La;


	// BC  u(x,0) = 0
	double dudy, f;
	for (int i = 0; i <= m_iMeshM + 1; i++)
	{
		m_pUij[i][0] = 0.;
	}

	// BC  u(x,1) = 0
	for (int i = 0; i <= m_iMeshM + 1; i++)
	{
		m_pUij[i][m_iMeshM + 1] = 0.;
	}

	// BC  u(0,y) = f(y)
	for (int j = 0; j <= m_iMeshM + 1; j++)
	{
		y = j * m_dH;
		m_pUij[0][j] = 0.;

	}

	// BC  u(1,y) = 0
	for (int j = 0; j <= m_iMeshM + 1; j++)
	{
		m_pUij[m_iMeshM + 1][j] = 0.;
	}

	// Inner grid
	for (int i = 1; i <= m_iMeshM; i++)
		for (int j = 1; j <= m_iMeshM; j++)
		{
			x = i * m_dH;
			y = j * m_dH;

			La = 1. / m_dH / m_dH * (m_pU0ij[i + 1][j] + m_pU0ij[i][j + 1] - 4.*m_pU0ij[i][j] + m_pU0ij[i - 1][j] + m_pU0ij[i][j - 1]);


			m_pUij[i][j] = 2.* m_pU0ij[i][j]- m_pU00ij[i][j] + m_dt * m_dt *La;;

		}


}


void Grid::SolveTwoDWave()
{
	m_dt = 1.e-4;  //m_dt must smaller than dx
	
	int iMax = (0.5 / m_dt) - 1;  //(1. / m_dt) - 1;
	//iMax = 1;

	TwoDWaveStep1();
	UpateTime();
	OutputResultAsVTK();

	for (int i = 1; i <= iMax; i++)
	{
		TwoDWaveStep2();
		UpateTime();
		//if ((i + 1) % 1000 == 0) OutputResultAsVTK();
	}

	OutputResultAsVTK();
}


void Grid::TwoDWaveStep1( )
{

	int iCount;
	double Error, MaxError;
	double x,y;
	double La;

	
		// BC  u(x,0) = 0
		double dudy, f;
		for (int i = 0; i <= m_iMeshM + 1; i++)
		{
			m_pUij[i][0] = 0.;
		}

		// BC  u(x,1) = 0
		for (int i = 0; i <= m_iMeshM + 1; i++)
		{
			m_pUij[i][m_iMeshM + 1] = 0.;
		}

		// BC  u(0,y) = f(y)
		for (int j = 0; j <= m_iMeshM + 1; j++)
		{
			y = j * m_dH;
			m_pUij[0][j] = 0.;

		}

		// BC  u(1,y) = 0
		for (int j = 0; j <= m_iMeshM + 1; j++)
		{
			m_pUij[m_iMeshM + 1][j] = 0.;
		}

		// Inner grid
		for (int i = 1; i <= m_iMeshM; i++)
			for (int j = 1; j <= m_iMeshM; j++)
			{
				x = i * m_dH;
				y = j * m_dH;

				La = 1./ m_dH/ m_dH*(m_pU0ij[i + 1][j] + m_pU0ij[i][j + 1] - 4.*m_pU0ij[i][j] + m_pU0ij[i - 1][j] + m_pU0ij[i][j - 1]);


				m_pUij[i][j] = m_pU0ij[i][j] + FunctionX(x)*FunctionX(y)*m_dt + m_dt * m_dt / 2.*La;;

			}

	

}

void Grid::UpateTime()
{
	m_iT = m_iT + 1;

	for (int i = 0; i <= m_iMeshM + 1; i++)
		for (int j = 0; j <= m_iMeshM + 1; j++)
		{
			{
				m_pU00ij[i][j] = m_pU0ij[i][j];
				m_pU0ij[i][j] = m_pUij[i][j];
			}
		}
}

void Grid::ExactSolution_Cooling()
{
	m_isExact = 1;
	//printf("Computing ExactSolution_Cooling \n");

	//m_TargetTime = 0.01;
	//m_Alpha = 1.;

	int max_series;
	max_series = 5;
	// all grid


	double x, y;
	double xx=0.;


	double m2_1 ,n2_1 ;
	double dTimeDecay;


	double MaxError = 1.e-30;
	double Error;



	//  ***  check m,n error
	
	m2_1 = 2.*max_series - 1.;
	n2_1 = 2.*max_series - 1.;
	dTimeDecay = exp(-(n2_1*n2_1 + m2_1 * m2_1)* m_Alpha * PI*PI * m_TargetTime);

	for (int i = 0; i <= m_iMeshM + 1; i++)
		for (int j = 0; j <= m_iMeshM + 1; j++)
		{
			{
				x = i * m_dH;
				y = j * m_dH;

				Error = fabs(1600. / m2_1 / n2_1 / PI / PI * sin(n2_1*PI*x) * sin(m2_1*PI*y) * dTimeDecay);
				if (Error > MaxError) MaxError = Error;
			}
		}

	//printf("m,n = %d, Time=%lg, Error= %lg\n", max_series, m_TargetTime, MaxError);

	
	for (int m = 1; m <= max_series; m++)
	{
		m2_1 = 2.*m - 1.;

		for (int n = 1; n <= max_series; n++)
		{
			n2_1 = 2.*n - 1.;

			dTimeDecay = exp(-(n2_1*n2_1 + m2_1 * m2_1)* m_Alpha * PI*PI * m_TargetTime);



			for (int i = 0; i <= m_iMeshM + 1; i++)
			{
				for (int j = 0; j <= m_iMeshM + 1; j++)
				{
					x = i * m_dH;
					y = j * m_dH;

					m_pUij[i][j] += 1600. / m2_1 / n2_1 / PI / PI * sin(n2_1*PI*x) * sin(m2_1*PI*y) * dTimeDecay;
				}
			}

		}

	}
	
	OutputResultAsVTK();

	//printf(" ExactSolution_Cooling Done \n");
}

void  Grid::SolveTriDiagonalMatrix(int iN, double *pA, double *pB, double *pC, double *pD)
{
	m_inTDM++;
	// 0~m+1  , iN = m+1


	pC[0] = pC[0] / pB[0];
	pD[0] = pD[0] / pB[0];

	for (int i = 1; i <= iN; i++)
	{
		pB[i] = pB[i] - pA[i] * pC[i - 1];
		pC[i] = pC[i] / pB[i];

		pD[i] = (pD[i] - pA[i] * pD[i - 1]) / pB[i];
	}

	for (int i = iN-1; i >= 0; i--)
	{
		pD[i] = pD[i] - pC[i] * pD[i + 1];
	}


}


void Grid::Set2DCoolingInitalCondiction()
{

	for (int i = 0; i <= m_iMeshM + 1; i++)
		for (int j = 0; j <= m_iMeshM + 1; j++)
		{
			{
				m_pUij[i][j]  = 100.;
				m_pU0ij[i][j] = 100.;
			
			}
		}

}

double Grid::ImplicitEulerLineX()
{
	double a_h2, One_dt;

	a_h2 = m_Alpha / m_dH / m_dH;
	One_dt = 1. / m_dt;

	double MaxError = 1.e-30;
	double Error;

	for (int j = 1; j <= m_iMeshM; j++)  // sweep y
	{
		for (int i = 1; i <= m_iMeshM; i++)  //only inner !
		{
			m_pA[i] = -a_h2;
			m_pB[i] = 4.*a_h2 + One_dt;
			m_pC[i] = -a_h2;
			m_pD[i] = m_pU0ij[i][j] * One_dt + a_h2 * (m_pUij[i][j - 1] + m_pUij[i][j + 1]) ;
		}

		m_pA[0] = 0.;        // BC
		m_pB[0] = 1.;
		m_pC[0] = 0.;
		m_pD[0] = 0.;

		m_pA[m_iMeshM + 1] = 0.;   // BC
		m_pB[m_iMeshM + 1] = 1.;
		m_pC[m_iMeshM + 1] = 0.;
		m_pD[m_iMeshM + 1] = 0.;

		SolveTriDiagonalMatrix(m_iMeshM + 1, m_pA, m_pB, m_pC, m_pD);

		for (int i = 0; i <= m_iMeshM + 1; i++)
		{

			Error = fabs(m_pUij[i][j] - m_pD[i]);    //Check error
			if (Error > MaxError) MaxError = Error;

			m_pUij[i][j] = m_pD[i];
		}

	}

	for (int i = 0; i <= m_iMeshM + 1; i++)   // BC
	{
		m_pUij[i][0] = 0.;
		m_pUij[i][m_iMeshM + 1] = 0.;
	}

	return MaxError;

}


double Grid::CrankNicolsonLineX()
{
	double a_h2, One_dt;

	a_h2 = m_Alpha / m_dH / m_dH;
	One_dt = 1. / m_dt;

	double MaxError = 1.e-30;
	double Error;

	for (int j = 1; j <= m_iMeshM; j++)  // sweep y
	{
		for (int i = 1; i <= m_iMeshM; i++)  //only inner !
		{
			m_pA[i] = -0.5*a_h2;
			m_pB[i] = 2.*a_h2 + One_dt;
			m_pC[i] = -0.5*a_h2;
			m_pD[i] = m_pU0ij[i][j] * One_dt + 0.5*a_h2 * (m_pUij[i][j - 1]  + m_pUij[i][j + 1]) + 0.5*a_h2 * (m_pU0ij[i + 1][j] + m_pU0ij[i][j + 1] - 4.*m_pU0ij[i][j] + m_pU0ij[i - 1][j] + m_pU0ij[i][j - 1]);
		}
	
		m_pA[0] = 0.;        // BC
		m_pB[0] = 1.;
		m_pC[0] = 0.;
		m_pD[0] = 0.;

		m_pA[m_iMeshM + 1] = 0.;   // BC
		m_pB[m_iMeshM + 1] = 1.;
		m_pC[m_iMeshM + 1] = 0.;
		m_pD[m_iMeshM + 1] = 0.;
	
		SolveTriDiagonalMatrix(m_iMeshM + 1, m_pA, m_pB, m_pC, m_pD);

		for (int i = 0; i <= m_iMeshM + 1; i++)
		{
			
			Error = fabs(m_pUij[i][j] - m_pD[i]);    //Check error
			if (Error > MaxError) MaxError = Error;

			m_pUij[i][j] = m_pD[i];
		}
	
	}

	for (int i = 0; i <= m_iMeshM + 1; i++)   // BC
	{
		m_pUij[i][0] = 0.;
		m_pUij[i][m_iMeshM + 1] = 0.;
	}

	return MaxError;

}

double Grid::CrankNicolsonLineY()
{
	double a_h2, One_dt;

	a_h2 = m_Alpha / m_dH / m_dH;
	One_dt = 1. / m_dt;

	double MaxError = 1.e-30;
	double Error;

	for (int i = 1; i <= m_iMeshM; i++)  // sweep x
	{
		for (int j = 1; j <= m_iMeshM; j++)  //only inner !
		{
			m_pA[j] = -0.5*a_h2;
			m_pB[j] = 2.*a_h2 + One_dt;
			m_pC[j] = -0.5*a_h2;
			m_pD[j] = m_pU0ij[i][j] * One_dt + 0.5*a_h2 * (m_pUij[i-1][j] + m_pUij[i+1][j]) + 0.5*a_h2 * (m_pU0ij[i + 1][j] + m_pU0ij[i][j + 1] - 4.*m_pU0ij[i][j] + m_pU0ij[i - 1][j] + m_pU0ij[i][j - 1]);
		}

		m_pA[0] = 0.;        // BC
		m_pB[0] = 1.;
		m_pC[0] = 0.;
		m_pD[0] = 0.;

		m_pA[m_iMeshM + 1] = 0.;   // BC
		m_pB[m_iMeshM + 1] = 1.;
		m_pC[m_iMeshM + 1] = 0.;
		m_pD[m_iMeshM + 1] = 0.;

		SolveTriDiagonalMatrix(m_iMeshM + 1, m_pA, m_pB, m_pC, m_pD);

		for (int j = 0; j <= m_iMeshM + 1; j++)
		{

			Error = fabs(m_pUij[i][j] - m_pD[j]);    //Check error
			if (Error > MaxError) MaxError = Error;

			m_pUij[i][j] = m_pD[j];
		}

	}

	for (int j = 0; j <= m_iMeshM + 1; j++)   // BC
	{
		m_pUij[0][j] = 0.;
		m_pUij[m_iMeshM + 1] [j]= 0.;
	}

	return MaxError;

}


void Grid::ADI_Step1()
{
	double a_h2, One_dt;

	a_h2 = m_Alpha / m_dH / m_dH;
	One_dt = 1. / m_dt;


	for (int j = 1; j <= m_iMeshM ; j++)  // sweep y
	{
		for (int i = 1; i <= m_iMeshM ; i++)  //only inner !
		{
			m_pA[i] = -0.5*a_h2;
			m_pB[i] = a_h2 + One_dt;
			m_pC[i] = -0.5*a_h2;
			m_pD[i] = m_pU0ij[i][j] * One_dt + 0.5*a_h2 * (m_pU0ij[i][j - 1] - 2.*m_pU0ij[i][j] + m_pU0ij[i][j + 1] );
		}
		//printf("1 \n");

		m_pA[0] = 0.;        // BC
		m_pB[0] = 1.;
		m_pC[0] = 0.;
		m_pD[0] = 0.;
		 
		m_pA[m_iMeshM + 1] = 0.;   // BC
		m_pB[m_iMeshM + 1] = 1.;
		m_pC[m_iMeshM + 1] = 0.;
		m_pD[m_iMeshM + 1] = 0.;

		//printf("2 \n");
		SolveTriDiagonalMatrix(m_iMeshM + 1, m_pA, m_pB, m_pC, m_pD);

		//printf("3 \n");

		for (int i = 0; i <= m_iMeshM + 1; i++)
		{
			m_pUij[i][j] = m_pD[i];
		}
		//printf("4 \n");
	}

	for (int i = 0; i <= m_iMeshM + 1; i++)   // BC
	{
		m_pUij[i][0] = 0.;
		m_pUij[i][m_iMeshM + 1] = 0.;
	}

	//printf("5 \n");

	    for (int i = 0; i <= m_iMeshM + 1; i++)
		for (int j = 0; j <= m_iMeshM + 1; j++)
		{
			{
				m_pU0ij[i][j] = m_pUij[i][j];
			}
		}

		/*printf("6 \n");
		getchar();*/
	
	//(m_pU0ij[i + 1][j] + m_pU0ij[i][j + 1] - 4.*m_pU0ij[i][j] + m_pU0ij[i - 1][j] + m_pU0ij[i][j - 1]);
}

void Grid::ADI_Step2()
{
	double a_h2, One_dt;

	a_h2 = m_Alpha / m_dH / m_dH;
	One_dt = 1. / m_dt;


	for (int i = 1; i <= m_iMeshM ; i++)  // sweep x
	{
		for (int j = 1; j <= m_iMeshM ; j++)  //only inner !
		{
			m_pA[j] = -0.5*a_h2;
			m_pB[j] = a_h2 + One_dt;
			m_pC[j] = -0.5*a_h2;
			m_pD[j] = m_pU0ij[i][j] * One_dt + 0.5*a_h2 * (m_pU0ij[i-1][j] - 2.*m_pU0ij[i][j] + m_pU0ij[i+1][j]);
		}


		m_pA[0] = 0.;  //BC
		m_pB[0] = 1.;
		m_pC[0] = 0.;
		m_pD[0] = 0.;

		m_pA[m_iMeshM + 1] = 0.;   //BC
		m_pB[m_iMeshM + 1] = 1.;
		m_pC[m_iMeshM + 1] = 0.;
		m_pD[m_iMeshM + 1] = 0.;


		SolveTriDiagonalMatrix(m_iMeshM + 1, m_pA, m_pB, m_pC, m_pD);

		for (int j = 0; j <= m_iMeshM + 1; j++)
		{
			m_pUij[i][j] = m_pD[j];
		}

	}

	for (int j = 0; j <= m_iMeshM + 1; j++)  //BC
	{
		m_pUij[0][j]            = 0.;
		m_pUij[m_iMeshM + 1][j] = 0.;
	}



	for (int i = 0; i <= m_iMeshM + 1; i++)
		for (int j = 0; j <= m_iMeshM + 1; j++)
		{
			{
				m_pU0ij[i][j] = m_pUij[i][j];
			}
		}

	/*printf("Step 2 -6 \n");
	getchar();*/


	//(m_pU0ij[i + 1][j] + m_pU0ij[i][j + 1] - 4.*m_pU0ij[i][j] + m_pU0ij[i - 1][j] + m_pU0ij[i][j - 1]);
}

void Grid::Solve2DCoolingByCrankNicolson(double dt)
{
	//printf("Solve2DCoolingByCrankNicolson \n");

	//m_dt = 1.e-6;  //m_dt must smaller than dx

	m_dt = dt;

	int iMax = (m_TargetTime / m_dt) ;//- 1
	int iMaxIter = 2000;

	double dError, dToler;

	dToler = 1.e-8;

	Set2DCoolingInitalCondiction();

	

	for (int iTime = 1; iTime <= iMax; iTime++)
	{
		for (int iter = 1; iter <= iMaxIter; iter++)
		{
			dError = CrankNicolsonLineX();
			//dError = CrankNicolsonLineY();
			
			if (dError <= dToler)
			{
				//printf("TimeLoop=%d, Iter=%d, Error =%lg ,LineByLine Solver \n", iTime, iter, dError);
				break;
			}		
			  
		} 
		
			for (int i = 0; i <= m_iMeshM + 1; i++)
			for (int j = 0; j <= m_iMeshM + 1; j++)
			{
				{
					m_pU0ij[i][j] = m_pUij[i][j];
				}
			}

		//if ((i + 1) % 1000 == 0) OutputResultAsVTK();
	}

	OutputResultAsVTK();

	//printf("Solve2DCoolingByCrankNicolson Done \n");

}

void Grid::Solve2DCoolingByImplicitEuler(double dt)
{
	//printf("Solve2DCoolingByCrankNicolson \n");

	//m_dt = 1.e-6;  //m_dt must smaller than dx

	m_dt = dt;

	int iMax = (m_TargetTime / m_dt);//- 1
	int iMaxIter = 2000;

	double dError, dToler;

	dToler = 1.e-8;

	Set2DCoolingInitalCondiction();



	for (int iTime = 1; iTime <= iMax; iTime++)
	{
		for (int iter = 1; iter <= iMaxIter; iter++)
		{
			dError = ImplicitEulerLineX();
			//dError = CrankNicolsonLineY();

			if (dError <= dToler)
			{
				//printf("TimeLoop=%d, Iter=%d, Error =%lg ,LineByLine Solver \n", iTime, iter, dError);
				break;
			}

		}

		for (int i = 0; i <= m_iMeshM + 1; i++)
			for (int j = 0; j <= m_iMeshM + 1; j++)
			{
				{
					m_pU0ij[i][j] = m_pUij[i][j];
				}
			}

		//if ((i + 1) % 1000 == 0) OutputResultAsVTK();
	}

	OutputResultAsVTK();

	//printf("Solve2DCoolingByCrankNicolson Done \n");

}

void Grid::Solve2DCoolingByADI(double dt)
{
	printf("Solve2DCoolingByADI \n");
	
	//m_dt = 1.e-4;  //m_dt must smaller than dx
	m_dt = dt;

	int iMax = (m_TargetTime / m_dt) ;


	Set2DCoolingInitalCondiction();

	//printf("Set2DCoolingInitalCondiction \n");

	for (int i = 1; i <= iMax; i++)
	{	
		ADI_Step1();

		ADI_Step2();

		//if ((i + 1) % 1000 == 0) OutputResultAsVTK();
	}

	OutputResultAsVTK();

	printf("Solve2DCoolingByADI Done \n");

}

void Grid::Solve2DCoolingByLOD(double dt)
{

	printf("Solve2DCoolingByLOD \n");


	//m_dt = 1.e-4;  //m_dt must smaller than dx

	m_dt = dt;

	int iMax = (m_TargetTime / m_dt) ; //- 1


	Set2DCoolingInitalCondiction();

	for (int i = 1; i <= iMax; i++)
	{
		LOD_Step1();

		LOD_Step2();

		//if ((i + 1) % 1000 == 0) OutputResultAsVTK();
	}

	OutputResultAsVTK();

	printf("Solve2DCoolingByLOD Done\n");

}

void Grid::LOD_Step1()
{
	double a_h2, One_dt;

	a_h2 = m_Alpha / m_dH / m_dH;
	One_dt = 1. / m_dt;


	for (int j = 0; j <= m_iMeshM + 1; j++)  // sweep y
	{
		for (int i = 1; i <= m_iMeshM ; i++)  //only inner !
		{
			m_pA[i] = -0.5*a_h2;
			m_pB[i] = a_h2 + One_dt;
			m_pC[i] = -0.5*a_h2;
			m_pD[i] = m_pU0ij[i][j] * One_dt + 0.5*a_h2 * (m_pU0ij[i-1][j] - 2.*m_pU0ij[i][j] + m_pU0ij[i+1][j]);
		}


		m_pA[0] = 0.;
		m_pB[0] = 1.;
		m_pC[0] = 0.;
		m_pD[0] = 0.;

		m_pA[m_iMeshM + 1] = 0.;
		m_pB[m_iMeshM + 1] = 1.;
		m_pC[m_iMeshM + 1] = 0.;
		m_pD[m_iMeshM + 1] = 0.;


		SolveTriDiagonalMatrix(m_iMeshM + 1, m_pA, m_pB, m_pC, m_pD);

		for (int i = 0; i <= m_iMeshM + 1; i++)
		{
			m_pUij[i][j] = m_pD[i];
		}

	}

	for (int i = 0; i <= m_iMeshM + 1; i++)
		for (int j = 0; j <= m_iMeshM + 1; j++)
		{
			{
				m_pU0ij[i][j] = m_pUij[i][j];
			}
		}


	//(m_pU0ij[i + 1][j] + m_pU0ij[i][j + 1] - 4.*m_pU0ij[i][j] + m_pU0ij[i - 1][j] + m_pU0ij[i][j - 1]);
}

void Grid::LOD_Step2()
{
	double a_h2, One_dt;

	a_h2 = m_Alpha / m_dH / m_dH;
	One_dt = 1. / m_dt;


	for (int i = 0; i <= m_iMeshM + 1; i++)  // sweep x
	{
		for (int j = 1; j <= m_iMeshM; j++)  //only inner !
		{
			m_pA[j] = -0.5*a_h2;
			m_pB[j] = a_h2 + One_dt;
			m_pC[j] = -0.5*a_h2;
			m_pD[j] = m_pU0ij[i][j] * One_dt + 0.5*a_h2 * (m_pU0ij[i][j-1] - 2.*m_pU0ij[i][j] + m_pU0ij[i][j+1]);
		}


		m_pA[0] = 0.;
		m_pB[0] = 1.;
		m_pC[0] = 0.;
		m_pD[0] = 0.;

		m_pA[m_iMeshM + 1] = 0.;
		m_pB[m_iMeshM + 1] = 1.;
		m_pC[m_iMeshM + 1] = 0.;
		m_pD[m_iMeshM + 1] = 0.;


		SolveTriDiagonalMatrix(m_iMeshM + 1, m_pA, m_pB, m_pC, m_pD);

		for (int j = 0; j <= m_iMeshM + 1; j++)
		{
			m_pUij[i][j] = m_pD[j];
		}

	}

	for (int i = 0; i <= m_iMeshM + 1; i++)
		for (int j = 0; j <= m_iMeshM + 1; j++)
		{
			{
				m_pU0ij[i][j] = m_pUij[i][j];
			}
		}


	//(m_pU0ij[i + 1][j] + m_pU0ij[i][j + 1] - 4.*m_pU0ij[i][j] + m_pU0ij[i - 1][j] + m_pU0ij[i][j - 1]);
}

void Grid::Solve1DInvisBurgers_nonConServative()
{
	printf("Solve1DInvisBurgers_nonConServative \n");

	double *pUj, *pU0j;

	pUj = m_pA;
	pU0j = m_pB;

	double dt_dx, x ;

	m_dt = 1.e-4;
	m_TargetTime = 1.;

	dt_dx = m_dt / m_dH;

	int iMax = (m_TargetTime / m_dt) ;


	char StrFileName[100] = { 0 };
	sprintf(StrFileName, "Solve1DInvisBurgers_nonConServative%d.csv", m_iMeshM);

	FILE *fp;
	fp = fopen(StrFileName, "w");


	fprintf(fp, "\n Time= %f \n", 0.);

	// Initial 
	Inital1DInvisBurgers();

	for (int j = 0; j <= m_iMeshM + 1; j++)
	{
		x = j * m_dH;
		fprintf(fp, "%f, %f \n", x, pUj[j]);
	}




	for (int iTime = 1; iTime <= iMax; iTime++)
	{


		for (int j = 1; j <= m_iMeshM; j++) // Inner
		{
			if (pU0j[j] > 0.)  pUj[j] = pU0j[j] - dt_dx * pU0j[j] * (pU0j[j] - pU0j[j - 1]);
			else              pUj[j] = pU0j[j] - dt_dx * pU0j[j] * (pU0j[j + 1] - pU0j[j]);

		}

		// periodic BC
		if (pU0j[0] > 0.)  pUj[0] = pU0j[0] - dt_dx * pU0j[0] * (pU0j[0] - pU0j[m_iMeshM]);  // -1 = m_iMeshM
		else              pUj[0] = pU0j[0] - dt_dx * pU0j[0] * (pU0j[1] - pU0j[0]);

		if (pU0j[m_iMeshM + 1] > 0.)  pUj[m_iMeshM + 1] = pU0j[m_iMeshM + 1] - dt_dx * pU0j[m_iMeshM + 1] * (pU0j[m_iMeshM + 1] - pU0j[m_iMeshM]);
		else                         pUj[m_iMeshM + 1] = pU0j[m_iMeshM + 1] - dt_dx * pU0j[m_iMeshM + 1] * (pU0j[1] - pU0j[m_iMeshM + 1]);  // m+2 = 1


		if ((iTime ) % 1000 == 0)
		{ 
			fprintf(fp, "\n Time= %f \n", iTime*m_dt);
			for (int j = 0; j <= m_iMeshM + 1; j++)
			{
				x = j * m_dH;
				fprintf(fp, "%f, %f \n", x, pUj[j]);
			}
		}

		for (int j = 0; j <= m_iMeshM + 1; j++)
		{
			pU0j[j] = pUj[j];
		}


	}


}

void Grid::Solve1DInvisBurgers_ConServative()
{
	printf("Solve1DInvisBurgers_ConServative \n");

	double *pUj, *pU0j;

	pUj = m_pA;
	pU0j = m_pB;

	double dt_dx, x;

	m_dt = 1.e-4;
	m_TargetTime = 1.;

	dt_dx = m_dt / m_dH;

	int iMax = (m_TargetTime / m_dt);


	char StrFileName[100] = { 0 };
	sprintf(StrFileName, "Solve1DInvisBurgers_ConServative%d.csv", m_iMeshM);

	FILE *fp;
	fp = fopen(StrFileName, "w");


	fprintf(fp, "\n Time= %f \n", 0.);

	// Initial 
	Inital1DInvisBurgers();

	for (int j = 0; j <= m_iMeshM + 1; j++)
	{
		x = j * m_dH;
		fprintf(fp, "%f, %f \n", x, pUj[j]);
	}




	for (int iTime = 1; iTime <= iMax; iTime++)
	{


		for (int j = 1; j <= m_iMeshM; j++) // Inner
		{
			if (pU0j[j] > 0.)  pUj[j] = pU0j[j] - dt_dx * 0.5 * (pU0j[j] * pU0j[j] - pU0j[j - 1] * pU0j[j - 1]);
			else              pUj[j] = pU0j[j] - dt_dx * 0.5 * (pU0j[j + 1] * pU0j[j + 1] - pU0j[j] * pU0j[j]);

		}

		// periodic BC
		if (pU0j[0] > 0.)  pUj[0] = pU0j[0] - dt_dx * 0.5 * (pU0j[0] * pU0j[0] - pU0j[m_iMeshM] * pU0j[m_iMeshM]);  // -1 = m_iMeshM
		else              pUj[0] = pU0j[0] - dt_dx * 0.5 * (pU0j[1] * pU0j[1] - pU0j[0] * pU0j[0]);

		if (pU0j[m_iMeshM + 1] > 0.)  pUj[m_iMeshM + 1] = pU0j[m_iMeshM + 1] - dt_dx * 0.5 * (pU0j[m_iMeshM + 1] * pU0j[m_iMeshM + 1] - pU0j[m_iMeshM] * pU0j[m_iMeshM]);
		else                         pUj[m_iMeshM + 1] = pU0j[m_iMeshM + 1] - dt_dx * 0.5 * (pU0j[1] * pU0j[1] - pU0j[m_iMeshM + 1] * pU0j[m_iMeshM + 1]);  // m+2 = 1


		if ((iTime) % 1000 == 0)
		{
			fprintf(fp, "\n Time= %f \n", iTime*m_dt);
			for (int j = 0; j <= m_iMeshM + 1; j++)
			{
				x = j * m_dH;
				fprintf(fp, "%f, %f \n", x, pUj[j]);
			}
		}

		for (int j = 0; j <= m_iMeshM + 1; j++)
		{
			pU0j[j] = pUj[j];
		}


	}


}

void Grid::Solve1DInvisBurgers_Godunov()
{
	printf("Solve1DInvisBurgers_Godunov \n");

	double *pUj, *pU0j;

	pUj = m_pA;
	pU0j = m_pB;

	double dt_dx, x;

	m_dt = 1.e-4;
	m_TargetTime = 10.;

	dt_dx = m_dt / m_dH;

	int iMax = (m_TargetTime / m_dt);


	char StrFileName[100] = { 0 };
	sprintf(StrFileName, "Solve1DInvisBurgers_Godunov%d.csv", m_iMeshM);

	FILE *fp;
	fp = fopen(StrFileName, "w");


	fprintf(fp, "\n Time= %f \n", 0.);

	// Initial 
	Inital1DInvisBurgers();

	for (int j = 0; j <= m_iMeshM + 1; j++)
	{
		x = j * m_dH;
		fprintf(fp, "%f, %f \n", x, pUj[j]);
	}




	for (int iTime = 1; iTime <= iMax; iTime++)
	{


		for (int j = 1; j <= m_iMeshM; j++) // Inner
		{
			 pUj[j] = pU0j[j] - dt_dx * ( Flux(j, j+1) - Flux(j-1, j) );
		}

		// periodic BC
		pUj[0] = pU0j[0] - dt_dx * ( Flux(0, 1) - Flux(m_iMeshM, 0)); // -1 = m_iMeshM

		pUj[m_iMeshM + 1] = pU0j[m_iMeshM + 1] - dt_dx * (Flux(m_iMeshM + 1, 1) - Flux(m_iMeshM, m_iMeshM + 1)  ); // m+2 = 1


		if ((iTime) % 1000 == 0)
		{
			fprintf(fp, "\n Time= %f \n", iTime*m_dt);
			for (int j = 0; j <= m_iMeshM + 1; j++)
			{
				x = j * m_dH;
				fprintf(fp, "%f, %f \n", x, pUj[j]);
			}
		}

		for (int j = 0; j <= m_iMeshM + 1; j++)
		{
			pU0j[j] = pUj[j];
		}


	}


}

void Grid::Inital1DInvisBurgers()
{
	double *pUj, *pU0j;

	pUj = m_pA;
	pU0j = m_pB;

	double x;

	for (int j = 0; j <= m_iMeshM + 1; j++)
	{
		x = j * m_dH;
		pUj[j] = pU0j[j] = 0.+ sin(2. * PI * x);  //  1.5

	}
}

double Grid::Flux(int iL, int iR)
{
	double dFlux=0.;
	double *pUj, *pU0j;

	pUj = m_pA;
	pU0j = m_pB;


	/*if ( pU0j[iL] > 0. && pU0j[iR] > 0.) return  0.5* pU0j[iL] * pU0j[iL];
	else if (pU0j[iL] < 0. && pU0j[iR] < 0.) return  0.5* pU0j[iR] * pU0j[iR];
	else return 0.;*/

	double s;

	s = 0.5*(pU0j[iL] + pU0j[iR]);

	if ( pU0j[iL] > 0. && s >0.) dFlux =   0.5* pU0j[iL] * pU0j[iL];
	if ( pU0j[iR] < 0. && s <0.)  dFlux =  0.5* pU0j[iR] * pU0j[iR];
	if ( pU0j[iL] < 0. && pU0j[iR] > 0.) dFlux = 0.;
	

	return dFlux;
}

void Grid::Solve2DLineByLine()
{
	for (int i = 0; i <= m_iMeshM + 1; i++)
	{
		m_pA[i] = 1.;
		m_pB[i] = -2.;
		m_pC[i] = 1.;
		m_pD[i] = 0.;
	}

	m_pA[0] = 0.;
	m_pB[0] = 1.;
	m_pC[0] = 0.;
	m_pD[0] = 0.;

	m_pA[m_iMeshM + 1] = 0.;
	m_pB[m_iMeshM + 1] = 1.;
	m_pC[m_iMeshM + 1] = 0.;
	m_pD[m_iMeshM + 1] = 1.;


	SolveTriDiagonalMatrix(m_iMeshM + 1, m_pA, m_pB, m_pC, m_pD);

	for (int i = 0; i <= m_iMeshM + 1; i++)
	{
		printf("i=%d, u=%lg \n", i, m_pD[i]);
	}

}

void Grid::SolveByJacobi(int iLoopMax, double dRelax)
{
	// Jacobi Loop
	//int iLoopMax = 1e6;

	m_iJacobiLoopMax = iLoopMax;
	double TOL = 1e-8;

	int iCount;
	double Error, MaxError;
	double y;

	for (iCount = 1; iCount <= m_iJacobiLoopMax; iCount++)
	{

		// BC  dudy(x,0) = 0
		double dudy, f;
		for (int i = 0; i <= m_iMeshM + 1; i++)
		{
			dudy = 0.;
			f = 0.;
			m_pUij[i][0] = m_pU0ij[i][1] - m_dH * dudy + m_dH * m_dH / 2.*f;
		}

		// BC  dudy(x,1) = 0
		for (int i = 0; i <= m_iMeshM + 1; i++)
		{
			dudy = 0.;
			f = 0.;
			m_pUij[i][m_iMeshM + 1] = m_pU0ij[i][m_iMeshM] - m_dH * dudy + m_dH * m_dH / 2.*f;
		}

		// BC  u(0,y) = f(y)
		for (int j = 0; j <= m_iMeshM + 1; j++)
		{
			y = j * m_dH;
			m_pUij[0][j] = cos(2.*PI*y);

			/*f = cos(2.*PI*y);
			if (f >= 0)
			{
				m_pUij[0][j] = 1;
			}
			else
			{
				m_pUij[0][j] = -1;
			}*/
			
		}

		// BC  u(1,y) = 0
		for (int j = 0; j <= m_iMeshM + 1; j++)
		{
			m_pUij[m_iMeshM + 1][j] = 0.;
		}

		// Inner grid
		for (int i = 1; i <= m_iMeshM; i++)
			for (int j = 1; j <= m_iMeshM; j++)
			{
				{
					m_pUij[i][j] = 0.25*(m_pU0ij[i + 1][j] + m_pU0ij[i][j + 1] + m_pU0ij[i - 1][j] + m_pU0ij[i][j - 1]);

					m_pUij[i][j] = m_pU0ij[i][j] + dRelax * (m_pUij[i][j] - m_pU0ij[i][j]);
				}
			}

		MaxError = 1.e-30;
		for (int i = 0; i <= m_iMeshM + 1; i++)
			for (int j = 0; j <= m_iMeshM + 1; j++)
			{
				{
					Error = fabs(m_pU0ij[i][j] - m_pUij[i][j]);
					if (Error > MaxError) MaxError = Error;

					m_pU0ij[i][j] = m_pUij[i][j];
				}
			}

		//printf("iCount=%d, MaxError=%lg \n", iCount, MaxError);

		if (MaxError < TOL) break;
	}

	m_dMaxError = MaxError;
	m_iIterations = iCount;

	if(iLoopMax >100) printf("    J MeshM=%d, Iterations=%d, MaxError=%lg \n", m_iMeshM, m_iIterations, m_dMaxError);
}

void Grid::OutputResultAsVTK()
{
	//output result
	double x, y;

	char StrFileName[100] = { 0 };

	//sprintf(StrFileName, "file_%d.vtk", m_iMeshM);

	if (m_isExact != 1)
	{
		sprintf(StrFileName, "file_FDM_Grid%d_Time%g.vtk", m_iMeshM, m_dt*m_iT);
	}
	else
	{
		sprintf(StrFileName, "file_Exact_Grid%d_Time%g.vtk", m_iMeshM, m_dt*m_iT);
	}
	

	FILE *fp;
	//fp = fopen("file.vtk", "w");
	fp = fopen(StrFileName, "w");
	fprintf(fp, "# vtk DataFile Version 3.0\n");
	fprintf(fp, "Result Field\n");
	fprintf(fp, "ASCII\n");
	fprintf(fp, "DATASET STRUCTURED_GRID\n");
	fprintf(fp, "DIMENSIONS %d %d %d\n", m_iMeshM + 2, m_iMeshM + 2, 1);
	fprintf(fp, "POINTS %d double\n", m_iNtotal);


	for (int i = 0; i <= m_iMeshM + 1; i++)
		for (int j = 0; j <= m_iMeshM + 1; j++)
		{
			{
				x = i * m_dH;
				y = j * m_dH;

				fprintf(fp, "%f %f %f\n", x, y, 0.);
			}
		}

	fprintf(fp, "\nPOINT_DATA %d\n", m_iNtotal);
	fprintf(fp, "SCALARS U double 1\n");
	fprintf(fp, "LOOKUP_TABLE default\n");

	for (int i = 0; i <= m_iMeshM + 1; i++)
		for (int j = 0; j <= m_iMeshM + 1; j++)
		{
			{
				fprintf(fp, "%f\n", m_pUij[i][j]);

			}
		}

	fclose(fp);//closing file.
}

void Grid::SolveByMultiGrid(int iLoopMax)
{
	// 2h grid
	int m_2h = (m_iMeshM - 1) / 2;

	double ** Rij_2h, ** R0ij_2h;
	Rij_2h = new double *[m_2h + 2];
	for (int i = 0; i <= m_2h + 1; i++) Rij_2h[i] = new double[m_2h + 2];

	R0ij_2h = new double *[m_2h + 2];
	for (int i = 0; i <= m_2h + 1; i++) R0ij_2h[i] = new double[m_2h + 2];

	for (int i = 0; i <= m_2h + 1; i++)
		for (int j = 0; j <= m_2h + 1; j++)
		{
			{
				Rij_2h[i][j] = 0.;
				R0ij_2h[i][j] = 0.;
			}
		}

	double TOL = 1e-8;

	int iCount;
	double Error, MaxError;
	
	for ( iCount = 1; iCount <= iLoopMax; iCount++)
	{

			for (int i = 0; i <= m_iMeshM + 1; i++)
				for (int j = 0; j <= m_iMeshM + 1; j++)
				{
					{

						m_pU00ij[i][j] = m_pU0ij[i][j];  //m_pU00ij used for error comparison
					}
				}


			//
			SolveByJacobi(3, 0.6666);

			// Compeute the residual of h grid


				for (int i = 1; i <= m_iMeshM; i++)     // only do inner
					for (int j = 1; j <= m_iMeshM; j++)
					{
						{
							m_pRij[i][j] = 0. - (m_pUij[i + 1][j] + m_pUij[i][j + 1] - 4.* (m_pUij[i][j]) + m_pUij[i - 1][j] + m_pUij[i][j - 1]);

							//printf("m_pRij[%d][%d]=%lg \n", i, j, m_pRij[i][j]);
						}
					}

				
				//for (int i = 0; i <= m_iMeshM + 1; i++) // BC  dudy(x,0) = 0
				//{
				//	m_pRij[i][0] = 0. - (-m_pUij[i][0] + m_pUij[i][1]);
				//}

				//
				//for (int i = 0; i <= m_iMeshM + 1; i++) // BC  dudy(x,1) = 0
				//{
				//	m_pRij[i][m_iMeshM + 1] = 0. - (-m_pUij[i][m_iMeshM + 1] + m_pUij[i][m_iMeshM]);
				//}

				// BC Dirichlet , residul =0


			// Coarsen the residual 
			for (int i = 0; i <= m_2h+1; i++)    
				for (int j = 0; j <= m_2h+1; j++)
				{
					{
						Rij_2h[i][j] = m_pRij[i * 2][j * 2];
						R0ij_2h[i][j] = Rij_2h[i][j];

						//printf("Rij_2h[%d][%d]=%lg \n", i, j, Rij_2h[i][j]);
					}
				}

			//Solve residual vector on 2h grid
			 for (int k = 1; k <= 3; k++)
			{

				for (int i = 1; i <= m_2h; i++)     // only do inner
					for (int j = 1; j <= m_2h; j++)
					{
						{
							Rij_2h[i][j] = 0.25*(R0ij_2h[i][j]  + R0ij_2h[i + 1][j] + R0ij_2h[i][j + 1] + R0ij_2h[i - 1][j] + R0ij_2h[i][j - 1]);
						}
					}
				
					//for (int i = 0; i <= m_2h + 1; i++) // BC  dudy(x,0) = 0
					//{
					//	Rij_2h[i][0] = (R0ij_2h[i][0] + R0ij_2h[i][1]);
					//}


					//for (int i = 0; i <= m_2h + 1; i++) // BC  dudy(x,1) = 0
					//{
					//	Rij_2h[i][m_2h + 1] = (R0ij_2h[i][m_2h + 1] + R0ij_2h[i][m_2h]);
					//}

				for (int i = 0; i <= m_2h+1; i++)    
					for (int j = 0; j <= m_2h+1; j++)
					{
						{

							R0ij_2h[i][j] = Rij_2h[i][j];

							//printf(" after Solve , Rij_2h[%d][%d]=%lg \n", i, j, Rij_2h[i][j]);
						}
					}
			}

			int iX, iY;
			// Interpolate back to h grid
			for (int i = 1; i <= m_iMeshM; i++)     // only do inner
				for (int j = 1; j <= m_iMeshM; j++)
				{

					if ((i % 2 == 0) && (j % 2 == 0))
					{
						iX = int(i / 2);
						iY = int(j / 2);
						m_pRij[i][j] = Rij_2h[iX][iY];
					}
					else if ((i % 2 == 0) && (j % 2 != 0))
					{
						iX = int(i / 2);
						iY = int(j / 2);
						m_pRij[i][j] = 0.5*(Rij_2h[iX][iY] + Rij_2h[iX][iY + 1]);
					}
					else if ((i % 2 != 0) && (j % 2 == 0))
					{
						iX = int(i / 2);
						iY = int(j / 2);
						m_pRij[i][j] = 0.5*(Rij_2h[iX][iY] + Rij_2h[iX + 1][iY]);
					}
					else
					{
						iX = int(i / 2);
						iY = int(j / 2);
						m_pRij[i][j] = 0.25*(Rij_2h[iX][iY] + Rij_2h[iX + 1][iY] + Rij_2h[iX][iY + 1] + Rij_2h[iX + 1][iY + 1]);
					}

					m_pUij[i][j] -= m_pRij[i][j];

					//printf(" Interpolate , m_pU0ij[%d][%d]=%lg \n", i, j, m_pU0ij[i][j]);
					//printf(" Interpolate , m_pUij[%d][%d]=%lg \n", i, j, m_pUij[i][j]);


					m_pU0ij[i][j] = m_pUij[i][j];
				}

			// iterations on h grid to smooth out error induced by interpolation
			SolveByJacobi(3, 0.6666);



			MaxError = 1.e-30;
			for (int i = 0; i <= m_iMeshM + 1; i++)
				for (int j = 0; j <= m_iMeshM + 1; j++)
				{
					{
						Error = fabs(m_pU00ij[i][j] - m_pUij[i][j]);
						if (Error > MaxError) MaxError = Error;
				
					}
				}

			printf("iCount=%d, MaxError=%lg \n", iCount, MaxError);

			if (MaxError < TOL) break;
	}

	m_dMaxError = MaxError;
	m_iIterations = iCount;

	//printf("MeshM=%d, Iterations=%d, MaxError=%lg \n", m_iMeshM, m_iIterations, m_dMaxError);
}

/*
void Grid::SolveByMultiGrid(int iLoopMax)
{
	// 2h grid
	int m_2h = (m_iMeshM - 1) / 2;

	double ** Rij_2h, ** R0ij_2h;
	Rij_2h = new double *[m_2h + 2];
	for (int i = 0; i <= m_iMeshM + 1; i++) Rij_2h[i] = new double[m_2h + 2];

	R0ij_2h = new double *[m_2h + 2];
	for (int i = 0; i <= m_iMeshM + 1; i++) R0ij_2h[i] = new double[m_2h + 2];

	for (int i = 0; i <= m_2h + 1; i++)
		for (int j = 0; j <= m_2h + 1; j++)
		{
			{
				Rij_2h[i][j] = 0.;
				R0ij_2h[i][j] = 0.;
			}
		}

	double TOL = 1e-8;

	int iCount;
	double Error, MaxError;

	for (iCount = 1; iCount <= iLoopMax; iCount++)
	{

		for (int i = 0; i <= m_iMeshM + 1; i++)
			for (int j = 0; j <= m_iMeshM + 1; j++)
			{
				{

					m_pU00ij[i][j] = m_pU0ij[i][j];  //m_pU00ij used for error comparison
				}
			}


		//
		SolveByJacobi(3, 0.6666);

		// Compeute the residual of h grid


		for (int i = 1; i <= m_iMeshM; i++)     // only do inner
			for (int j = 1; j <= m_iMeshM; j++)
			{
				{
					m_pRij[i][j] = 0. - (m_pUij[i + 1][j] + m_pUij[i][j + 1] - 4.* (m_pUij[i][j]) + m_pUij[i - 1][j] + m_pUij[i][j - 1]);

					//printf("m_pRij[%d][%d]=%lg \n", i, j, m_pRij[i][j]);
				}
			}



		// Coarsen the residual 
		for (int i = 1; i <= m_2h; i++)     // only do inner
			for (int j = 1; j <= m_2h; j++)
			{
				{
					Rij_2h[i][j] = m_pRij[i * 2][j * 2];
					R0ij_2h[i][j] = Rij_2h[i][j];

					//printf("Rij_2h[%d][%d]=%lg \n", i, j, Rij_2h[i][j]);
				}
			}

		//Solve residual vector on 2h grid
		for (int k = 1; k <= 3; k++)
		{

			for (int i = 1; i <= m_2h; i++)     // only do inner
				for (int j = 1; j <= m_2h; j++)
				{
					{
						Rij_2h[i][j] = 0.25*(R0ij_2h[i][j] + R0ij_2h[i + 1][j] + R0ij_2h[i][j + 1] + R0ij_2h[i - 1][j] + R0ij_2h[i][j - 1]);
					}
				}

			for (int i = 1; i <= m_2h; i++)     // only do inner
				for (int j = 1; j <= m_2h; j++)
				{
					{

						R0ij_2h[i][j] = Rij_2h[i][j];

						//printf(" after Solve , Rij_2h[%d][%d]=%lg \n", i, j, Rij_2h[i][j]);
					}
				}
		}

		int iX, iY;
		// Interpolate back to h grid
		for (int i = 1; i <= m_iMeshM; i++)     // only do inner
			for (int j = 1; j <= m_iMeshM; j++)
			{

				if ((i % 2 == 0) && (j % 2 == 0))
				{
					iX = int(i / 2);
					iY = int(j / 2);
					m_pRij[i][j] = Rij_2h[iX][iY];
				}
				else if ((i % 2 == 0) && (j % 2 != 0))
				{
					iX = int(i / 2);
					iY = int(j / 2);
					m_pRij[i][j] = 0.5*(Rij_2h[iX][iY] + Rij_2h[iX][iY + 1]);
				}
				else if ((i % 2 != 0) && (j % 2 == 0))
				{
					iX = int(i / 2);
					iY = int(j / 2);
					m_pRij[i][j] = 0.5*(Rij_2h[iX][iY] + Rij_2h[iX + 1][iY]);
				}
				else
				{
					iX = int(i / 2);
					iY = int(j / 2);
					m_pRij[i][j] = 0.25*(Rij_2h[iX][iY] + Rij_2h[iX + 1][iY] + Rij_2h[iX][iY + 1] + Rij_2h[iX + 1][iY + 1]);
				}

				m_pUij[i][j] -= m_pRij[i][j];

				//printf(" Interpolate , m_pU0ij[%d][%d]=%lg \n", i, j, m_pU0ij[i][j]);
				//printf(" Interpolate , m_pUij[%d][%d]=%lg \n", i, j, m_pUij[i][j]);


				m_pU0ij[i][j] = m_pUij[i][j];
			}

		// iterations on h grid to smooth out error induced by interpolation
//		SolveByJacobi(3, 0.6666);



		MaxError = 1.e-30;
		for (int i = 0; i <= m_iMeshM + 1; i++)
			for (int j = 0; j <= m_iMeshM + 1; j++)
			{
				{
					Error = fabs(m_pU00ij[i][j] - m_pUij[i][j]);
					if (Error > MaxError) MaxError = Error;

				}
			}

		printf("iCount=%d, MaxError=%lg \n", iCount, MaxError);

		if (MaxError < TOL) break;
	}

	m_dMaxError = MaxError;
	m_iIterations = iCount;

	//printf("MeshM=%d, Iterations=%d, MaxError=%lg \n", m_iMeshM, m_iIterations, m_dMaxError);
}

*/

void Grid::NewSpectralMemory()
{

	m_pDNij = new double *[m_iChebyN + 1];
	for (int i = 0; i <= m_iChebyN; i++) m_pDNij[i] = new double[m_iChebyN + 1];

	m_pDN2ij = new double *[m_iChebyN + 1];
	for (int i = 0; i <= m_iChebyN; i++) m_pDN2ij[i] = new double[m_iChebyN + 1];


	for (int i = 0; i <= m_iChebyN; i++)
		for (int j = 0; j <= m_iChebyN; j++)
		{
			{
				m_pDNij[i][j] = 0.;
				m_pDN2ij[i][j] = 0.;
			}
		}

	m_pDN2ij_Tu = new double *[m_iChebyN - 1];
	for (int i = 0; i <= m_iChebyN - 2; i++) m_pDN2ij_Tu[i] = new double[m_iChebyN - 1];

	m_pIij_Tu = new double *[m_iChebyN - 1];
	for (int i = 0; i <= m_iChebyN - 2; i++) m_pIij_Tu[i] = new double[m_iChebyN - 1];

	for (int i = 0; i <= m_iChebyN - 2; i++)
		for (int j = 0; j <= m_iChebyN - 2; j++)
		{
			{
				m_pDN2ij_Tu[i][j] = 0.;
				m_pIij_Tu[i][j] = 0.;

				if (i == j) m_pIij_Tu[i][j] = 1.;
			}
		}

	m_iChebyInner2DGrids = (m_iChebyN - 1)*(m_iChebyN - 1);

	m_pLNij = new double *[m_iChebyInner2DGrids];
	for (int i = 0; i < m_iChebyInner2DGrids; i++) m_pLNij[i] = new double[m_iChebyInner2DGrids];

	for (int i = 0; i < m_iChebyInner2DGrids; i++)
		for (int j = 0; j < m_iChebyInner2DGrids; j++)
		{
			{
				m_pLNij[i][j] = 0.;
			}
		}


	m_pUInner = new double[m_iChebyInner2DGrids];

	m_pUtInner = new double[m_iChebyInner2DGrids];

	m_pLaU = new double[m_iChebyInner2DGrids];

	m_pLaUt = new double[m_iChebyInner2DGrids];

	m_pLaULaU = new double[m_iChebyInner2DGrids];

	for (int i = 0; i < m_iChebyInner2DGrids; i++)
	{
		m_pUInner[i] = 0.;
		m_pLaU[i] = 0.;

		m_pLaULaU[i] = 0.;

		m_pUtInner[i] = 0.;

		m_pLaUt[i] = 0.;
	}
}
double Grid::ChebyGridX(int j)
{

	return cos(j * PI / m_iChebyN);
}

void Grid::SetChebyshevDifferMatrix()
{
	double xi, xj;

	for (int i = 1; i <= m_iChebyN-1; i++)
		for (int j = 1; j <= m_iChebyN-1; j++)
		{
			{
				xi = ChebyGridX(i);
				xj = ChebyGridX(j);

				if (i == j)
				{
					m_pDNij[i][j] = -xj / (2.*(1 - xj * xj));
				}
				else
				{
					m_pDNij[i][j] = pow(-1, i + j) / (xi - xj);
				}
			}
		}

	for (int i = 1; i <= m_iChebyN - 1; i++) 
	{
		xi = ChebyGridX(i);
		
		m_pDNij[i][0] = -1. / 2. *   pow(-1, i) / (1 - xi);  // left 

		m_pDNij[i][m_iChebyN] = 1. / 2. *   pow(-1, i+ m_iChebyN) / (1 + xi);  // right 

	}

	for (int j = 1; j <= m_iChebyN - 1; j++)
	{
		xj = ChebyGridX(j);

		m_pDNij[0][j] = 2. *   pow(-1, j) / (1 - xj);  // Top 

		m_pDNij[m_iChebyN][j] = -2. *   pow(-1, j + m_iChebyN) / (1 + xj);  // buttom 

	}

	m_pDNij[0][0] = (2.*m_iChebyN*m_iChebyN + 1.) / 6.;
	m_pDNij[m_iChebyN][m_iChebyN] = -(2.*m_iChebyN*m_iChebyN + 1.) / 6.;

	m_pDNij[0][m_iChebyN] = 1. / 2. *  pow(-1, m_iChebyN);
	m_pDNij[m_iChebyN][0] = -1. / 2. *  pow(-1, m_iChebyN);


	////Debug printf
	//printf("\n");
	//for (int i = 0; i <= m_iChebyN; i++)
	//{	
	//	for (int j = 0; j <= m_iChebyN; j++)
	//		{

	//			printf("%lg ", m_pDNij[i][j]);
	//			//m_pDNij[i][j] = 0.;		
	//		}
	//		printf("\n");
	//}	
}

void Grid::SetChebyshevDNSquare()
{
	for (int i = 0; i <= m_iChebyN; i++)
		{
			for (int j = 0; j <= m_iChebyN; j++)
			{
				for (int k = 0; k <= m_iChebyN; k++)
				{
					m_pDN2ij[i][j] += m_pDNij[i][k] * m_pDNij[k][j];
				}
			}
		}

	//// Debug print
	//printf("\n");
	//for (int i = 0; i <= m_iChebyN ; i++)

	//{
	//	for (int j = 0; j <= m_iChebyN ; j++)
	//	{

	//		printf("%lg ", m_pDN2ij[i][j]);
	//	}
	//	printf("\n");
	//}


	for (int i = 0; i <= m_iChebyN - 2; i++)
		for (int j = 0; j <= m_iChebyN - 2; j++)
		{
			{
				m_pDN2ij_Tu[i][j] = m_pDN2ij[i+1][j+1];
				
			}
		}

	//// Debug print
	//printf("\n");
	//for (int i = 0; i <= m_iChebyN - 2; i++)
	//	
	//	{
	//	for (int j = 0; j <= m_iChebyN - 2; j++)
	//		{
	//			
	//			printf("%lg ", m_pDN2ij_Tu[i][j]);
	//		}
	//		printf("\n");
	//	}

}

void Grid::SetChebyLaplacian()
{
	int rowA, colA;
	int rowB, colB;

	int startRow, startCol;

	rowA = colA = m_iChebyN - 1;
	rowB = colB = m_iChebyN - 1;

	for (int i = 0; i < rowA; i++) {
		for (int j = 0; j < colA; j++) {
			startRow = i * rowB;
			startCol = j * colB;
			for (int k = 0; k < rowB; k++) {
				for (int l = 0; l < colB; l++) {
					m_pLNij[startRow + k][startCol + l] = m_pIij_Tu[i][j] * m_pDN2ij_Tu[k][l]  + m_pDN2ij_Tu[i][j] * m_pIij_Tu[k][l];
				}
			}
		}
	}

	//// Debug print
	//printf("\n");
	//for (int i = 0; i < m_iChebyInner2DGrids; i++)
	//	{
	//		for (int j = 0; j < m_iChebyInner2DGrids; j++)
	//		{
	//			
	//			printf("%lg ", m_pLNij[i][j]);
	//		}
	//		printf("\n");
	//	}

}


void Grid::SolveTwoDWavebySpectral()
{
	m_dt = 2.5e-4;  //

	int iMax = (0.5 / m_dt) - 1; //(1. / m_dt) - 1
	//iMax = 1;

	NewSpectralMemory();

	 SetChebyshevDifferMatrix();
	 SetChebyshevDNSquare();
	 SetChebyLaplacian();

	 ComputeChebyLaplacianUt();

	 UpdateUInnerForSpectral();
	 TwoDWaveStep1bySpectral();
	 UpateTime();
	 OutputResultAsVTKSpectral();


	for (int i = 1; i <= iMax; i++)
	{
		UpdateUInnerForSpectral();
		TwoDWaveStep2bySpectral();
		UpateTime();
		//if ((i + 1) % 1000 == 0) OutputResultAsVTK();
	}

	OutputResultAsVTKSpectral();

}

void Grid::UpdateUInnerForSpectral()
{
	int k;
	// Inner grid
	for (int i = 1; i <= m_iMeshM; i++)
	{
		for (int j = 1; j <= m_iMeshM; j++)
		{
			k = (i - 1) + (j - 1)*(m_iChebyN - 1);

			m_pUInner[k] = m_pU0ij[i][j];

		}

	}
			
}

void Grid::ComputeChebyLaplacianUt()
{

	int k;
	// Inner grid
	for (int i = 1; i <= m_iMeshM; i++)
	{
		for (int j = 1; j <= m_iMeshM; j++)
		{
			k = (i - 1) + (j - 1)*(m_iChebyN - 1);

			m_pUtInner[k] = FunctionXSpectral(ChebyGridX(i))*FunctionXSpectral(ChebyGridX(j));

		}

	}


	double xx = 0.;

	for (int i = 0; i < m_iChebyInner2DGrids; i++)
	{
		xx = 0.;
		for (int j = 0; j < m_iChebyInner2DGrids; j++)
		{
			xx += m_pLNij[i][j] * m_pUtInner[j];
		}
		m_pLaUt[i] = xx;
	}
}

void Grid::ComputeChebyLaplacianU()
{
	double xx = 0.;

	for (int i = 0; i < m_iChebyInner2DGrids; i++)
	{
		xx = 0.;
		for (int j = 0; j < m_iChebyInner2DGrids; j++)
		{
			xx += m_pLNij[i][j] * m_pUInner[j];
		}
		m_pLaU[i] = xx;
	}

	for (int i = 0; i < m_iChebyInner2DGrids; i++)
	{
		xx = 0.;
		for (int j = 0; j < m_iChebyInner2DGrids; j++)
		{
			xx += m_pLNij[i][j] * m_pLaU[j];
		}
		m_pLaULaU[i] = xx;
	}
}

void Grid::TwoDWaveStep1bySpectral() 
{

	ComputeChebyLaplacianU();


	int iCount;
	double Error, MaxError;
	
	int k;

	// BC  u(x,0) = 0
	
	for (int i = 0; i <= m_iMeshM + 1; i++)
	{
		m_pUij[i][0] = 0.;
	}

	// BC  u(x,1) = 0
	for (int i = 0; i <= m_iMeshM + 1; i++)
	{
		m_pUij[i][m_iMeshM + 1] = 0.;
	}

	// BC  u(0,y) = f(y)
	for (int j = 0; j <= m_iMeshM + 1; j++)
	{
		m_pUij[0][j] = 0.;

	}

	// BC  u(1,y) = 0
	for (int j = 0; j <= m_iMeshM + 1; j++)
	{
		m_pUij[m_iMeshM + 1][j] = 0.;
	}

	// Inner grid
	for (int i = 1; i <= m_iMeshM; i++)
	{
		for (int j = 1; j <= m_iMeshM; j++)
		{
		
			k = (i - 1) + (j - 1)*(m_iChebyN - 1);

			m_pUij[i][j] = m_pU0ij[i][j] + FunctionXSpectral(ChebyGridX(i))*FunctionXSpectral(ChebyGridX(j))*m_dt 
				+ 1./12. *m_dt * m_dt *m_pLaUt[k]
				+ m_dt * m_dt / 2.*4.*(m_pLaU[k] - m_dt * m_dt /12.*m_pLaULaU[k]);

		}

	}
		
}

void Grid::TwoDWaveStep2bySpectral()
{
	ComputeChebyLaplacianU();


	int iCount;
	double Error, MaxError;

	int k;

	// BC  u(x,0) = 0

	for (int i = 0; i <= m_iMeshM + 1; i++)
	{
		m_pUij[i][0] = 0.;
	}

	// BC  u(x,1) = 0
	for (int i = 0; i <= m_iMeshM + 1; i++)
	{
		m_pUij[i][m_iMeshM + 1] = 0.;
	}

	// BC  u(0,y) = f(y)
	for (int j = 0; j <= m_iMeshM + 1; j++)
	{
		m_pUij[0][j] = 0.;

	}

	// BC  u(1,y) = 0
	for (int j = 0; j <= m_iMeshM + 1; j++)
	{
		m_pUij[m_iMeshM + 1][j] = 0.;
	}

	// Inner grid
	for (int i = 1; i <= m_iMeshM; i++)
	{
		for (int j = 1; j <= m_iMeshM; j++)
		{

			k = (i - 1) + (j - 1)*(m_iChebyN - 1);

			m_pUij[i][j] = 2.* m_pU0ij[i][j] - m_pU00ij[i][j] + m_dt * m_dt *4.*(m_pLaU[k] - m_dt * m_dt / 12.*m_pLaULaU[k]);
		}

	}


}

void Grid::OutputResultAsVTKSpectral()
{
	//output result
	double x, y;

	char StrFileName[100] = { 0 };

	//sprintf(StrFileName, "file_%d.vtk", m_iMeshM);

	sprintf(StrFileName, "file_Grid%d_Time%g.vtk", m_iMeshM, m_dt*m_iT);

	FILE *fp;
	//fp = fopen("file.vtk", "w");
	fp = fopen(StrFileName, "w");
	fprintf(fp, "# vtk DataFile Version 3.0\n");
	fprintf(fp, "Result Field\n");
	fprintf(fp, "ASCII\n");
	fprintf(fp, "DATASET STRUCTURED_GRID\n");
	fprintf(fp, "DIMENSIONS %d %d %d\n", m_iMeshM + 2, m_iMeshM + 2, 1);
	fprintf(fp, "POINTS %d double\n", m_iNtotal);


	for (int i = 0; i <= m_iMeshM + 1; i++)
		for (int j = 0; j <= m_iMeshM + 1; j++)
		{
			{
				x = ChebyGridX(i);  //i * m_dH;
				y = ChebyGridX(j); //j * m_dH;

				fprintf(fp, "%f %f %f\n", x, y, 0.);
			}
		}

	fprintf(fp, "\nPOINT_DATA %d\n", m_iNtotal);
	fprintf(fp, "SCALARS U double 1\n");
	fprintf(fp, "LOOKUP_TABLE default\n");

	for (int i = 0; i <= m_iMeshM + 1; i++)
		for (int j = 0; j <= m_iMeshM + 1; j++)
		{
			{
				fprintf(fp, "%f\n", m_pUij[i][j]);

			}
		}

	fclose(fp);//closing file.
}