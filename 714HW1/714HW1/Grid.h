#pragma once
class Grid
{
public:
	Grid();
	Grid(int m);
	~Grid();

	void SolveByJacobi(int iLoopMax, double dRelax); // HW1
	void SolveByMultiGrid(int iLoopMax);             // HW2
    void OutputResultAsVTK();                        // HW1


	void Solve1DInvisBurgers_nonConServative();  //HW4
	void Solve1DInvisBurgers_ConServative();
	void Solve1DInvisBurgers_Godunov();
	double Flux(int iL, int iR);
	void Inital1DInvisBurgers();

	void ExactSolution_Cooling();  //Final Project
	void SolveTriDiagonalMatrix(int iN, double *pA, double *pB, double *pC, double *pD);
	void Solve2DLineByLine();
	void Set2DCoolingInitalCondiction();
	void Solve2DCoolingByADI(double dt);
	void Solve2DCoolingByLOD(double dt);
	void Solve2DCoolingByCrankNicolson(double dt);
	void Solve2DCoolingByImplicitEuler(double dt);

	void ADI_Step1();
	void ADI_Step2();

	void LOD_Step1();
	void LOD_Step2();

	double ImplicitEulerLineX();
	double CrankNicolsonLineX();
	double CrankNicolsonLineY();

    
	void SolveTwoDWave();
	void TwoDWaveStep1();
	void TwoDWaveStep2();
	void UpateTime();

	//Spectral Method
	int m_iChebyN;
	int m_iChebyInner2DGrids;
	double *m_pUInner, *m_pLaU, *m_pLaULaU;
	double *m_pUtInner, *m_pLaUt;
	double **m_pDNij, **m_pDN2ij;           //0~N , dimension N+1
	double **m_pDN2ij_Tu, **m_pIij_Tu;          //0~N-2, dimension N-1   
	double **m_pLNij;                           // 0~ (N-1)^2-1,  dimension (N-1)^2   
	double ChebyGridX(int j);
	void NewSpectralMemory();
	void SetChebyshevDifferMatrix();
	void SetChebyshevDNSquare();
	void SetChebyLaplacian();
	void UpdateUInnerForSpectral();
	void SolveTwoDWavebySpectral();
	void ComputeChebyLaplacianU();
	void ComputeChebyLaplacianUt();
	void TwoDWaveStep1bySpectral();
	void TwoDWaveStep2bySpectral();
	void OutputResultAsVTKSpectral();

	int m_iMeshM;
	int m_iNtotal;
	int m_iT;
	int m_iIterations;
	int m_iJacobiLoopMax;
	double m_dH;
	double m_dMaxError;
	double m_dt;
	double m_TargetTime;

	double m_Alpha;  //Final Project
	int m_inTDM;
	int m_isExact;

	double **m_pUij, **m_pU0ij , **m_pRij;
	double **m_pU00ij;

	double *m_pA, *m_pB, *m_pC, *m_pD; // for TriDiagonalMatrix

};

