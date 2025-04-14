/*
Matrix Cx(5);
	Matrix Cy(5);
	Matrix Cz(5);
	Cx(1)=1;Cx(2)=2;Cx(3)=3;Cx(4)=4;Cx(5)=5;
	Cy(1)=-1;Cy(2)=-2;Cy(3)=-3;Cy(4)=-4;Cy(5)=-5;
	Cz(1)=0;Cz(2)=1;Cz(3)=2;Cz(4)=3;Cz(5)=4;
	double t = 2;
	double N = 5;
	double Ta = 1;
	double Tb = 3;
	
	Matrix A(3);
    A = Cheb3D(t,N,Ta,Tb,Cx,Cy,Cz);
*/