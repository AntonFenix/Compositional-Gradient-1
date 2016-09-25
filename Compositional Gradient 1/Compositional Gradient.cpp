#include <iostream>
using namespace std;

#define M_PI (3.141592653589793)
#define M_2PI (2.*M_PI)

class Compositional {	
	
	public:
		const int n = 15;
		const double g = 9.8;
		int m;
		bool cond = true;
		double p = 31.6;
		double* M;
		double* w;
		double* pc;
		double* Tc;
		double* z;	
		double Z;
		double** coeff;
		double *CPVAP_A, *CPVAP_B, *CPVAP_C, *CPVAP_D;

		Compositional() {
			M = new double[n];
			w = new double[n];
			pc = new double[n];
			Tc = new double[n];
			z = new double[n];

			coeff = new double*[n];
			for (int i = 0; i < n; i++) {
				coeff[i] = new double[n];
			}

			A = new double[n];
			a = new double[n];
			B = new double[n];
			b = new double[n];
			fug = new double[n];
			psi = new double[n];
		}
		
		void constantEOS();
		double* cordano(double *x);
		double* fugacity();
		double** J();
		double* F(double *y);

	private: 		
		const double err = 0.0001;
		const double h1 = -3000;
		const double h2 = -2800;
		double T = 334;
		double R = 8.31*1000;
		double Am, am, *A, *a;
		double Bm, bm, *B, *b;
		double *fug;	
		double *psi;

};

void Compositional::constantEOS() {
	double Tr;
	double m;
	double alpha;

	am = 0;
	bm = 0;
	for (int i = 0; i < n; i++) {
		Tr = T / Tc[i];

		if (w[i] > 0.4) {
			m = 0.3796 + 1.485*w[i] - 0.1644*pow(w[i], 2) + 0.01667*pow(w[i], 3);
		}
		else {
			m = 0.37464 + 1.54226*w[i] - 0.26992*pow(w[i], 2);
		}

		alpha = pow(1 + m*(1 - sqrt(Tr)), 2);
		a[i] = (0.45724*pow(R*Tc[i], 2)*alpha) / pc[i];
		//cout << "a[" << i + 1 << "] = " << a[i] << endl;
		b[i] = (0.07780*R*Tc[i]) / pc[i];
	}

	for (int i = 0; i < n; i++) {
		//cout << "b[" << i + 1 << "] = " << b[i] << endl;
		bm += z[i] * b[i];
		for (int j = 0; j < n; j++) {
			am += z[i] * z[j] * (1 - coeff[i][j]) * sqrt(a[i] * a[j]);
		}
	}
	Am = (am*p) / pow(R*T, 2);
	Bm = (bm*p) / (R*T);
	/*cout << "am = " << am << endl;
	cout << "bm = " << bm << endl;
	cout << "Am = " << Am << endl;
	cout << "Bm = " << Bm << endl << endl;*/
}

double* Compositional::cordano(double *x) {
	double b, c, d;
	b = -(1 - Bm);
	c = Am - 3 * pow(Bm, 2) - 2 * Bm;
	d = -(Am*Bm - pow(Bm, 2) - pow(Bm, 3));

	double p, q, D, F;

	p = c - pow(b, 2) / 3.;
	q = 2. / 27. * pow(b, 3) - 1. / 3. * b*c + d;
	D = pow((p / 3), 3) + pow((q / 2), 2);

	if (D < 0) {
		if (q < 0) {
			F = atan(sqrt(-D) / (-q / 2));
		}
		else if (q > 0) {
			F = atan(sqrt(-D) / (-q / 2)) + M_PI;
		}
		else {
			F = M_PI / 2;
		}

		x[0] = 2 * sqrt(-p / 3)*cos(F / 3) - b / 3;
		x[1] = 2 * sqrt(-p / 3)*cos(F / 3 + M_2PI / 3) - b / 3;
		x[2] = 2 * sqrt(-p / 3)*cos(F / 3 + 2.*M_2PI / 3) - b / 3;
		cout << "x[0] = " << x[0] << endl;
		cout << "x[1] = " << x[1] << endl;
		cout << "x[2] = " << x[2] << endl << endl;
	}
	else if (D > 0) {

		if ((-q / 2 - sqrt(D)) < 0)
			x[0] = pow((-q / 2 + sqrt(D)), 1. / 3.) - pow((q / 2 + sqrt(D)), 1. / 3.) - b / 3;
		else
			x[0] = pow((-q / 2 + sqrt(D)), 1. / 3.) + pow((-q / 2 - sqrt(D)), 1. / 3.) - b / 3;
		x[1] = x[0] - b / 3;
		x[2] = x[0] - b / 3;
		cout << "x[0] = " << x[0] << endl;
		cout << "x[1] = " << x[1] << endl;
		cout << "x[2] = " << x[2] << endl << endl;
	}
	else {
		x[0] = 2 * pow(-q / 2, 1. / 3.) - b / 3;
		x[1] = pow(-q / 2, 1. / 3.) - b / 3;
		x[2] = pow(-q / 2, 1. / 3.) - b / 3;
		cout << "x[0] = " << x[0] << endl;
		cout << "x[1] = " << x[1] << endl;
		cout << "x[2] = " << x[2] << endl << endl;
	}

	return x;
}

double* Compositional::fugacity()
{
	
	double C0, C1, C2, C3, C4, C5;
	for (int i = 0; i < n; i++) {
		C3 = 0;
		for (int j = 0; j < n; j++) {
			C3 += z[j] * (1 - coeff[i][j]) * sqrt(a[i] * a[j]);
		}

		C0 = log(z[i] * p);
		C1 = log(Z - Bm);
		C2 = (b[i] / bm)*(Z - 1);
		C4 = (Am / (Bm * 2 * sqrt(2)))*(b[i] / bm - (2 / am)*C3);
		C5 = log((Z + (1 + sqrt(2))*Bm) / (Z + (1 - sqrt(2))*Bm));
		fug[i] = C0 - C1 + C2 + C4*C5; 
	}

	if (m == 0) {
		for (int i = 0; i < n; i++) {
			psi[i] = fug[i] + M[i] * g * (h1 - h2) / (R*T);
			cout << "fug[" << i << "] = " << fug[i] << endl;
			cout << "Doveska[" << i << "] = " << M[i] * g * (h1 - h2) / (R*T) << endl;
		}
	}
	return fug;
}

double** Compositional::J() {
	double Cm = (1 + sqrt(2))*Bm, Dm = (1 - sqrt(2))*Bm;
	double cm = (1 + sqrt(2))*bm, dm = (1 - sqrt(2))*bm;

	double* d_ln_f_dp = new double[n], dF_dZ, dF_dz, dF_dp, dZ_dp, dZ_dz;
	double delta, Ml, Ql;
	double C1, C2;
	double Q1, Q2, Q3, Q4;
	double Z1, Z2, Z3;
	double f1, f2, f3;
	double P1, P2, P3, P4, P5, P6;

	double *C = new double[n], *D = new double[n], *B = new double[n], *c = new double[n], *d = new double[n];
	for (int i = 0; i < 15; i++) {
		B[i] = b[i] * p / (R*T);
		C[i] = (1 + sqrt(2))*B[i];
		D[i] = (1 - sqrt(2))*B[i];
		d[i] = (1 - sqrt(2))*b[i];
		c[i] = (1 + sqrt(2))*b[i];

	}

	double **d_ln_f_dz = new double*[n];
	for (int i = 0; i < n; i++) {
		d_ln_f_dz[i] = new double[n];
	}

	dF_dZ = 3 * pow(Z, 2) + 2 * Z*(Cm + Dm - Bm - 1) + (Am - Bm*Cm + Cm*Dm - Bm*Dm - Dm - Cm);
	dF_dp = (1 / p)*((Cm + Dm - Bm)*pow(Z, 2) + (Am - 2 * Bm*Cm + 2 * Cm*Dm - 2 * Bm*Dm - Dm - Cm)*Z - (3 * Bm*Cm*Dm + 2 * Cm*Dm + 2 * Am*Bm));
	dZ_dp = -dF_dp / dF_dZ;

	for (int i = 0; i < n; i++) {

		//Производная логарифма летучести по давлению
		P6 = 0;
		for (int j = 0; j < n; j++) {
			P6 = z[j] * (1 - coeff[i][j]) * sqrt(a[i] * a[j]);
		}

		P1 = (dZ_dp - bm / (R*T))*(1 + B[i] / (Z - Bm)) - b[i] / (R*T);
		P2 = 2 * P6 / am - (c[i] - d[i]) / (cm - dm);//+++
		P3 = (dZ_dp + cm / (R*T)) / (Z + Cm) - (dZ_dp + dm / (R*T)) / (Z + Dm);
		P4 = (c[i] / (R*T))*(Z + Cm) - C[i] * (dZ_dp + cm / (R*T));
		P5 = (d[i] / (R*T))*(Z + Dm) - D[i] * (dZ_dp + dm / (R*T));

		d_ln_f_dp[i] = 1 / p - (1 / (Z - Bm))*P1 - (Am / (Cm - Dm))*(P2*P3 + (1 / pow((Z + Cm), 2))*P4 - (1 / pow((Z + Dm), 2))*P5);
		//cout << "d_ln_f_dp[" << i << "] = " << d_ln_f_dp[i] << endl;
		

		//Производная логарифма летучести по компонентному составу
		for (int j = 0; j < n; j++) {
			if (i == j) delta = 1;
			else		delta = 0;

			C1 = 0; C2 = 0;
			for (int k = 0; k < 15; k++) {
				C1 += z[k] * (1 - coeff[i][k]) * sqrt(a[i] * a[k]);
				C2 += z[k] * (1 - coeff[j][k]) * sqrt(a[j] * a[k]);
			}

			Z1 = c[j] + d[j] - b[j];
			Z2 = (2 / (R*T))*C2 - Cm*b[j] - Bm*c[j] + Dm*c[j] + Cm*d[j] - Dm*b[j] - Bm*d[j] - d[j] - c[j];
			Z3 = -Cm*Dm*b[j] - Bm*Dm*c[j] - Bm*Cm*d[j] - Dm*c[j] - Cm*d[j] - Am*b[j] - (2 * Bm / (R*T))*C2;
			dF_dz = (p / (R*T))*(Z1*pow(Z, 2) + Z2*Z + Z3);

			dZ_dz = -dF_dz / dF_dZ;

			Ml = (1 / (Cm - Dm))*((2 * p / pow(R*T, 2))*C2 - (Am / (Cm - Dm))*(C[j] - D[j]));

			Q1 = 2 * ((1 - coeff[i][j]) * sqrt(a[i] * a[j])) / am - 4 * C1*C2 / pow(am, 2) + (c[i] - d[i])*(c[j] - d[j]) / pow((cm - dm), 2);
			Q2 = 2 * C1 / am - (c[i] - d[i]) / (cm - dm);
			Q3 = (1 / (Z + Cm))*(dZ_dz + C[j]) - (1 / (Z + Dm))*(dZ_dz + D[j]);
			Q4 = -(C[i] / pow((Z + Cm), 2))*(dZ_dz + C[j]) + (D[i] / pow((Z + Dm), 2))*(dZ_dz + D[j]);

			Ql = Q1*log((Z + Cm) / (Z + Dm)) + Q2*Q3 + Q4;

			f1 = (1 / (Z - Bm))*(dZ_dz - B[j])*(1 + B[i] / (Z - Bm));
			f2 = 2 * C1 / am - (c[i] - d[i]) / (cm - dm);
			f3 = C[i] / (Z + Cm) - D[i] / (Z + Dm);

			d_ln_f_dz[i][j] = delta / z[j] - f1 - Ml*(f2*log((Z + Cm) / (Z + Dm)) + f3) - Am*Ql / (Cm - Dm);
			//cout << "d_ln_f_dz[" << i << "][" << j << "] = " << d_ln_f_dz[i][j] << endl;

		}
	}

	double **F = new double*[n];
	for (int i = 0; i < n; i++) {
		F[i] = new double[n];
	}

	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			if (j == 0) F[i][j] = d_ln_f_dp[i];
			else	    F[i][j] = d_ln_f_dz[i][j] - d_ln_f_dz[i][0];
		}
	}
	delete C, D, B, c, d, d_ln_f_dp, d_ln_f_dz;

	return F;
}

double* Compositional::F(double* y) {
	
	for (int i = 0; i < n; i++) {
		y[i] = -(fug[i] - psi[i]);
		if (abs(y[i]) > err) cond = true;
	}
	return y;
}

// Вывод системы уравнений
void sysout(double **a, double *y, int n) {
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++){
			cout << a[i][j] << "*x" << j;
			if (j < n - 1) {
				cout << " + ";
			}
		}
		cout << " = " << y[i] << endl;
	}
	return;
}
double * gauss(double **a, double *y, int n) {
	double *x, max;
	int k, index;
	const double eps = 0.00001;  // точность
	x = new double[n];
	k = 0;
	while (k < n) {
		// Поиск строки с максимальным a[i][k]
		max = abs(a[k][k]);
		index = k;
		for (int i = k + 1; i < n; i++) {
			if (abs(a[i][k]) > max) {
				max = abs(a[i][k]);
				index = i;
			}
		}
		// Перестановка строк
		if (max < eps) {
			// нет ненулевых диагональных элементов
			cout << "Решение получить невозможно из-за нулевого столбца ";
			cout << index << " матрицы A" << endl;
			return 0;
		}
		for (int j = 0; j < n; j++) {
			double temp = a[k][j];
			a[k][j] = a[index][j];
			a[index][j] = temp;
		}
		double temp = y[k];
		y[k] = y[index];
		y[index] = temp;
		// Нормализация уравнений
		for (int i = k; i < n; i++) {
			double temp = a[i][k];
			if (abs(temp) < eps) continue; // для нулевого коэффициента пропустить
			for (int j = 0; j < n; j++) {
				a[i][j] = a[i][j] / temp;
			}
			y[i] = y[i] / temp;
			if (i == k)  continue; // уравнение не вычитать само из себя
			for (int j = 0; j < n; j++) {
				a[i][j] = a[i][j] - a[k][j];
			}
			y[i] = y[i] - y[k];
		}
		k++;
	}
	// обратная подстановка
	for (k = n - 1; k >= 0; k--) {
		x[k] = y[k];
		for (int i = 0; i < k; i++) {
			y[i] = y[i] - a[i][k] * x[k];
		}
	}
	return x;
}

void main() {
	system("chcp 1251");
	system("mode con cols=160 lines=5000");

	double M[15] = { 28.02, 44.01, 16.04, 30.07, 44.09, 58.12, 58.12, 72.15, 72.15, 85.01883, 110.4571, 157.9535, 231.2177, 338.1838, 500. };
	double w[15] = { 0.0450, 0.2310, 0.0115, 0.0908, 0.1454, 0.1756, 0.1928, 0.2273, 0.2510, 0.2922825, 0.2581414, 0.3432671, 0.8143111, 1.001544, 1.258908 };
	double pc[15] = { 3.399, 7.382, 4.604, 4.880, 4.249, 3.648, 3.797, 3.381, 3.369, 3.06011, 2.505358, 1.941391, 1.382302, 0.9461165, 0.6445802 };
	double Tc[15] = { 126.3, 304.2, 190.6, 305.4, 369.8, 408.2, 425.2, 460.4, 469.7, 504.8923, 570.4422, 646.6256, 732.3722, 818.5594, 904.8372 };
	double z[15] = { 0.05009623, 0.000399622, 0.7484794, 0.08561797, 0.03473826, 0.00445777, 0.007988606, 0.001606675, 0.001410077, 0.02081841, 0.02668487, 0.01381346, 0.003451272, 0.000403741, 0.0000336818 };
	double coeff[15][15] = { { 0., 0., 0.025, 0.01, 0.09, 0.095, 0.095, 0.1, 0.11, 0.11, 0.11, 0.11, 0.11, 0.11, 0.11 }, { 0., 0., 0.105, 0.13, 0.125, 0.12, 0.115, 0.115, 0.115, 0.115, 0.115, 0.115, 0.115, 0.115, 0.115 }, { 0.025, 0.105, 0., 0.002994556, 0.008343353, 0.01402962, 0.01402962, 0.01952127, 0.01952127, 0.02425953, 0.06545012, 0.06480181, 0.05242532, 0.09938945, 0.179825 }, { 0.01, 0.13, 0.002994556, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0. }, { 0.09, 0.125, 0.008343353, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0. }, { 0.095, 0.12, 0.01402962, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0. }, { 0.095, 0.115, 0.01402962, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0. }, { 0.1, 0.115, 0.01952127, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0. }, { 0.11, 0.115, 0.01952127, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0. }, { 0.11, 0.115, 0.02425953, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0. }, { 0.11, 0.115, 0.06545012, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0. }, { 0.11, 0.115, 0.06480181, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0. }, { 0.11, 0.115, 0.05242532, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0. }, { 0.11, 0.115, 0.09938945, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0. }, { 0.11, 0.115, 0.179825, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0. } };
	double CPVAP_A[9] = {3.115E+1, 1.980E+1, 1.925E+1, 5.409E+0, -4.224E+0, -1.390E+0, 9.487E+0, -9.525E+0, -3.626E+0};
	double CPVAP_B[9] = {-1.357E-2, 7.344E-2, 5.213E-2, 1.781E-1, 3.063E-1, 3.847E-1, 3.313E-1, 5.066E-1, 4.873E-1};
	double CPVAP_C[9] = {2.680E-5, -5.602E-5, 1.197E-5, -6.938E-5, -1.586E-4, -1.846E-4, -1.108E-4, -2.729E-4, -2.580E-4};
	double CPVAP_D[9] = {-1.168E-8, 1.715E-8, -1.132E-8, 8.713E-9, 3.215E-8, 2.895E-8, -2.822E-9, 5.723E-8, -2.822E-9};
	
	Compositional compos;
	compos.m = 0;

	compos.M = M;
	compos.w = w;
	compos.pc = pc;
	compos.Tc = Tc;
	compos.z = z;
	compos.CPVAP_A = CPVAP_A;
	compos.CPVAP_B = CPVAP_B;
	compos.CPVAP_C = CPVAP_C;
	compos.CPVAP_D = CPVAP_D;
	for (int i = 0; i < compos.n; i++) {
		for (int j = 0; j < compos.n; j++) {
			compos.coeff[i][j] = coeff[i][j];
		}
	}
	
	double** J, *F, *S;
	double sum = 0;
	double diff = 0;
	double *y = new double[compos.n], *x = new double[3];

	do{
		compos.cond = false;

		compos.constantEOS();

		compos.Z = compos.cordano(x)[0];
		cout << "Z = " << compos.Z << endl << endl;

		compos.fugacity(); 

		J = compos.J(); 
		cout << endl; 
		F = compos.F(y); 

		cout << endl;
		S = gauss(J, F, compos.n);
		
		cout << endl;

		sum = 0;
		compos.p += S[0];
		
		cout << "p = " << compos.p << endl << endl;
		for (int i = 1; i < compos.n; i++) {
			compos.z[i] += S[i];
			sum += compos.z[i];
		}
		compos.z[0] = 1 - sum;

		for (int i = 0; i < compos.n; i++){
			cout << "z[" << i << "]=" << compos.z[i] << endl;
		}

		compos.m++;
		cout << "m = " << compos.m << endl;

		cout << "//////////////////////////////////////////" << endl;
		cout << "//////////////////////////////////////////" << endl << endl;
				
	} while (compos.cond);
}