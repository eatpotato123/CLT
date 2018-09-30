#pragma omp parallel for \
		default(shared) \
		firstprivate(t) \
		private(B_par,b_par,curl_b_par,grad_B_par,B_value_par,psi_par,\
				deltaB_par,deltaE_par,deltaPhi_par,deltaA_par_par,grad_deltaA_par_par,grad_deltaPhi_par,\
				deltaE_X_par,deltaE_Y_par,\
				Va,RHS_V,Vold,\
                Ua,RHS_U,Uold,lambda,\
				Wa,RHS_W,Wold,\
				Xa,RHS_X,Xold,\
				timestep,coeff,told,\
				OUT_OF_BOUDARY,\
				weight_line,weight_square,weight_line_aux,weight_line_auxiliary,\
				i,j,R,Z,phi,B_s,B_ss,W_per,W_par,mu,\
				E,v,v_d,\
				grad_p_phi,dXdt1,dv_pardt1,\
				p_phi,\
				df0dt,dp_phidt,df_0dP_phi,dEdt,df_0dE,\
				deltaE_local,deltaE_X_local,deltaE_Y_local,deltaB_local,deltaPhi_local,deltaA_par_local,grad_deltaA_par_local,grad_deltaPhi_local,grad_deltaPhi_X_local,grad_deltaPhi_Y_local,\
                grad_deltaPhi_X_par,grad_deltaPhi_Y_par,\
                dZ_of_step1,Z_reverse_sign,\
                COS_NPHI_PLUS_OMEGA_T,SIN_NPHI_PLUS_OMEGA_T)
for (int p = 0; p<N_p; p++)
{

	OUT_OF_BOUDARY = false;
	mu = marker[p].mu;
	for (int RK4 = 1; RK4 <= 4; RK4++)
	{
		switch (RK4) {
		case 1:
		{
			Xold = marker[p].X;
			Vold = marker[p].v_par;
			Wold = marker[p].w;
			Xa = marker[p].X;
			Va = marker[p].v_par;
			Wa = marker[p].w;
			if (SLOWING_DOWN)
			{
				Uold = marker[p].v_per;
				Ua = marker[p].v_per;
			}
			told = t;
			timestep = 0.5*dt;
			coeff = 1.0 / 6.0;
		}break;
		case 2:
		{
			timestep = 0.5*dt;
			coeff = 1.0 / 3.0;
		}break;
		case 3:
		{
			timestep = dt;
			coeff = 1.0 / 3.0;
		}break;
		case 4:
		{
			timestep = dt;
			coeff = 1.0 / 6.0;
		}break;
		}

		i = (int)((marker[p].X[0] - R0 + a*a_b) / dR);
		j = (int)((marker[p].X[1] + a*a_b) / dZ);



		R = R0 - a*a_b + i*dR;
		Z = -a*a_b + j*dZ;



		weight_line_auxiliary[0][0] = (marker[p].X[0] - R) / dR;
		weight_line_auxiliary[1][0] = (marker[p].X[1] - Z) / dZ;

		weight_line_auxiliary[0][1] = 1.0 - weight_line_auxiliary[0][0];
		weight_line_auxiliary[1][1] = 1.0 - weight_line_auxiliary[1][0];


		for (int ii = 0; ii<2; ii++)
			for (int jj = 0; jj<2; jj++)
				weight_square[ii][jj] = weight_line_auxiliary[0][1 - ii] * weight_line_auxiliary[1][1 - jj];


		B_par = 0.0;
		curl_b_par = 0.0;
		grad_B_par = 0.0;
		B_value_par = 0.0;

		for (int ii = 0; ii<2; ii++)
			for (int jj = 0; jj<2; jj++)
			{
				B_value_par += B_value[i + ii][j + jj] * weight_square[ii][jj];
				for (int d = 0; d<2; d++)
				{
					B_par[d] += B[d][i + ii][j + jj] * weight_square[ii][jj];
					grad_B_par[d] += grad_B[d][i + ii][j + jj] * weight_square[ii][jj];
				}

				for (int d = 0; d<3; d++)
					curl_b_par[d] += curl_b[d][i + ii][j + jj] * weight_square[ii][jj];
			}

		if (NUMERICAL_EQUILIBRIUM)
		{
			double g_eq_par = 0.0;
			for (int ii = 0; ii<2; ii++)
				for (int jj = 0; jj<2; jj++)
					g_eq_par += g_eq[i + ii][j + jj] * weight_square[ii][jj];

			B_par[2] = g_eq_par / marker[p].X[0];
		}
		else
			B_par[2] = B0*R0 / marker[p].X[0];

		grad_B_par[2] = 0.0;


		for (int d = 0; d<3; d++)
			b_par[d] = B_par[d] / B_value_par;

		grad_deltaA_par_par = 0.0;
		grad_deltaPhi_par = 0.0;
		deltaA_par_par = 0.0;
		for (int s = 0; s<NUM_MODE; s++)
		{
			COS_NPHI_PLUS_OMEGA_T = cos(n[s] * marker[p].X[2] - omega_A[s] * t);
			SIN_NPHI_PLUS_OMEGA_T = sin(n[s] * marker[p].X[2] - omega_A[s] * t);

			for (int d = 0; d<3; d++)
				for (int I = i; I <= i + 1; I++)
					for (int J = j; J <= j + 1; J++)
					{
						if (d != 2)
							grad_deltaA_par_local[d][I - i][J - j] = X[s] * grad_deltaA_par_cos[s][d][I][J] * COS_NPHI_PLUS_OMEGA_T - X[s] * grad_deltaA_par_sin[s][d][I][J] * SIN_NPHI_PLUS_OMEGA_T + Y[s] * grad_deltaA_par_sin[s][d][I][J] * COS_NPHI_PLUS_OMEGA_T + Y[s] * grad_deltaA_par_cos[s][d][I][J] * SIN_NPHI_PLUS_OMEGA_T;
						else
							grad_deltaA_par_local[d][I - i][J - j] = (-X[s] * grad_deltaA_par_cos[s][d][I][J] * SIN_NPHI_PLUS_OMEGA_T - X[s] * grad_deltaA_par_sin[s][d][I][J] * COS_NPHI_PLUS_OMEGA_T - Y[s] * grad_deltaA_par_sin[s][d][I][J] * SIN_NPHI_PLUS_OMEGA_T + Y[s] * grad_deltaA_par_cos[s][d][I][J] * COS_NPHI_PLUS_OMEGA_T)*n[s] / marker[p].X[0];
						if (d != 2)
							grad_deltaPhi_local[d][I - i][J - j] = X[s] * grad_deltaPhi_cos[s][d][I][J] * COS_NPHI_PLUS_OMEGA_T - X[s] * grad_deltaPhi_sin[s][d][I][J] * SIN_NPHI_PLUS_OMEGA_T + Y[s] * grad_deltaPhi_sin[s][d][I][J] * COS_NPHI_PLUS_OMEGA_T + Y[s] * grad_deltaPhi_cos[s][d][I][J] * SIN_NPHI_PLUS_OMEGA_T;
						else
							grad_deltaPhi_local[d][I - i][J - j] = (-X[s] * grad_deltaPhi_cos[s][d][I][J] * SIN_NPHI_PLUS_OMEGA_T - X[s] * grad_deltaPhi_sin[s][d][I][J] * COS_NPHI_PLUS_OMEGA_T - Y[s] * grad_deltaPhi_sin[s][d][I][J] * SIN_NPHI_PLUS_OMEGA_T + Y[s] * grad_deltaPhi_cos[s][d][I][J] * COS_NPHI_PLUS_OMEGA_T)*n[s] / marker[p].X[0];

					}
			for (int I = i; I <= i + 1; I++)
				for (int J = j; J <= j + 1; J++)
					deltaA_par_local[I - i][J - j] = X[s] * deltaA_par_cos[s][I][J] * COS_NPHI_PLUS_OMEGA_T - X[s] * deltaA_par_sin[s][I][J] * SIN_NPHI_PLUS_OMEGA_T + Y[s] * deltaA_par_sin[s][I][J] * COS_NPHI_PLUS_OMEGA_T + Y[s] * deltaA_par_cos[s][I][J] * SIN_NPHI_PLUS_OMEGA_T;




			for (int ii = 0; ii<2; ii++)
				for (int jj = 0; jj<2; jj++)
				{
					for (int d = 0; d<3; d++)
					{
						grad_deltaA_par_par[d] += grad_deltaA_par_local[d][ii][jj] * weight_square[ii][jj];
						grad_deltaPhi_par[d] += grad_deltaPhi_local[d][ii][jj] * weight_square[ii][jj];
					}
					deltaA_par_par += deltaA_par_local[ii][jj] * weight_square[ii][jj];
				}


		}
		deltaB_par = times(grad_deltaA_par_par, b_par) + deltaA_par_par*curl_b_par;
		deltaE_par = (b_par*grad_deltaPhi_par)*b_par - grad_deltaPhi_par;

		B_s = B_par + deltaB_par + alpha*m*marker[p].v_par / e*curl_b_par;
		B_ss = B_s*b_par;



		if (!FULL_F && !RES && !COF)
		{
			/*----------------------------------------------------------------------|
			|                                                                       |
			|               delta-f  weight equation : dw/dt = - 1/g*df_0dt         |
			|                                                                       |
			-----------------------------------------------------------------------*/

			R = marker[p].X[0];
			Z = marker[p].X[1];

			E = 0.5*m*pow(marker[p].v_par, 2.0) + mu*B_value_par;
			v = sqrt(2 * E / species[0].mass);

			dXdt1 = 1.0 / B_ss*(marker[p].v_par*deltaB_par - times(b_par, deltaE_par));
			dv_pardt1 = e / m / B_ss*(B_s*deltaE_par / alpha - mu / e*deltaB_par*grad_B_par);

			v_d = 1.0 / (e*B_ss)*(alpha*m*pow(marker[p].v_par, 2.0)*curl_b_par + mu*times(b_par, grad_B_par));

			if (SHAFRANOV)
			{
				psi_par = 0;
				for (int ii = 0; ii<2; ii++)
					for (int jj = 0; jj<2; jj++)
						psi_par += psi[i + ii][j + jj] * weight_square[ii][jj];
			}
			else
			{
				double r2 = (R - R0)*(R - R0) + Z*Z;
				psi_par = psi_0*(q_0 + 1) + 0.5*sqrt(4 * pow(psi_0*(q_0 + 1), 2.0) - 4 * psi_0*B0*(r2 - a*a));
			}


			if (NUMERICAL_EQUILIBRIUM)
			{
				particle_vector<double> grad_g_eq_par;
				grad_g_eq_par = 0.0;
				for (int ii = 0; ii<2; ii++)
					for (int jj = 0; jj<2; jj++)
						for (int d = 0; d<3; d++)
							grad_g_eq_par[d] += grad_g_eq[d][ii][jj] * weight_square[ii][jj];

				grad_p_phi[0] = m*marker[p].v_par / B_value_par*grad_g_eq_par[0] - m*marker[p].v_par*R*B_par[2] / B_value_par / B_value_par*grad_B_par[0] - e*B_par[1] * R / alpha;
				grad_p_phi[1] = m*marker[p].v_par / B_value_par*grad_g_eq_par[1] - m*marker[p].v_par*R*B_par[2] / B_value_par / B_value_par*grad_B_par[1] + e*B_par[0] * R / alpha;
				grad_p_phi[2] = 0.0;
			}
			else
			{
				grad_p_phi[0] = -m*marker[p].v_par*R*B_par[2] / B_value_par / B_value_par*grad_B_par[0] - e*B_par[1] * R / alpha;
				grad_p_phi[1] = -m*marker[p].v_par*R*B_par[2] / B_value_par / B_value_par*grad_B_par[1] + e*B_par[0] * R / alpha;
				grad_p_phi[2] = 0.0;
			}



			dp_phidt = dXdt1*grad_p_phi + dv_pardt1*m*R*B_par[2] / B_value_par;
			dEdt = e*v_d*deltaE_par + e*marker[p].v_par*deltaB_par*deltaE_par / B_ss;


			/*----------------------------------------------------------------------------------|
			|                                                                                   |
			|                           dw/dt = -1/g*df_0/dt                                    |
			|                                                                                   |
			-----------------------------------------------------------------------------------*/

			p_phi = m*marker[p].v_par*marker[p].X[0] * (B_par[2] / B_value_par) + e*psi_par / alpha;

			if (AVERAGE_PSI)
			{
				double bracket_psi;
				int sgn_v_parallel = marker[p].v_par>0 ? 1 : -1;
				if ((1 - mu*B0 / E) > 0)
					bracket_psi = p_phi / species[0].charge - species[0].mass / species[0].charge*sgn_v_parallel*v*R0*sqrt(1 - mu*B0 / E);
				else
					bracket_psi = p_phi / species[0].charge;

				double f0 = c_f*1.0 / (pow(v, 3) + v_c*v_c*v_c)*(1 + erf((v_0 - v) / (0.2*v_A)))*exp(-(bracket_psi / (c_1*Delta_psi)));

				df_0dP_phi = -1.0 / (e*c_1*Delta_psi)*f0;

				if ((1 - mu*B0 / E) > 0)
					df_0dE = -3 / m*sqrt(2 * E / m) / (pow(2 * E / m, 1.5) + pow(v_c, 3.0))*f0
					- 1 / (sqrt(2 * PI*E*m)*0.1*v_A)*exp(-(pow((v_0 - v) / (0.2*v_A), 2.0)))*c_f*1.0 / (pow(v, 3) + v_c*v_c*v_c)*exp(-(bracket_psi / (c_1*Delta_psi)))
					+ sgn_v_parallel*R0 / (species[0].charge*c_1*Delta_psi*sqrt(v*v - 2 * mu*B0 / species[0].mass))*f0;
				else
					df_0dE = -3 / m*sqrt(2 * E / m) / (pow(2 * E / m, 1.5) + pow(v_c, 3.0))*f0
					- 1 / (sqrt(2 * PI*E*m)*0.1*v_A)*exp(-(pow((v_0 - v) / (0.2*v_A), 2.0)))*c_f*1.0 / (pow(v, 3) + v_c*v_c*v_c)*exp(-(bracket_psi / (c_1*Delta_psi)));
			}
			else
			{
				df_0dP_phi = -1.0 / (e*c_1*Delta_psi)*c_f*1.0 / (pow(v, 3) + v_c*v_c*v_c)*(1 + erf((v_0 - v) / (0.2*v_A)))*exp(-(p_phi / (e*c_1*Delta_psi)));

				df_0dE = -3 / m*sqrt(2 * E / m) / (pow(2 * E / m, 1.5) + pow(v_c, 3.0))*c_f*1.0 / (pow(v, 3) + v_c*v_c*v_c)*(1 + erf((v_0 - v) / (0.2*v_A)))*exp(-(p_phi / (e*c_1*Delta_psi)))
					- 1 / (sqrt(2 * PI*E*m)*0.1*v_A)*exp(-(pow((v_0 - v) / (0.2*v_A), 2.0)))*c_f*1.0 / (pow(v, 3) + v_c*v_c*v_c)*exp(-(p_phi / (e*c_1*Delta_psi)));
			}

			df0dt = dp_phidt*df_0dP_phi + dEdt*df_0dE;


			if (abs(df_0dE) > 1e10)
			{
				cout << "df_0dE : " << df_0dE << " v : " << v << " (1e-20 + 1+erf((v_0-v)/(0.2*v_A))) : " << (1e-20 + 1 + erf((v_0 - v) / (0.2*v_A))) << " exp(-(pow((v_0-v)/(0.2*v_A),2.0))) : " << exp(-(pow((v_0 - v) / (0.2*v_A), 2.0))) << " myid : " << myid << endl;
				exit(0);
			}

		}


		if (FULL_F || RES || COF)
		{
			if (SHAFRANOV)
			{
				psi_par = 0;
				for (int ii = 0; ii<2; ii++)
					for (int jj = 0; jj<2; jj++)
						psi_par += psi[i + ii][j + jj] * weight_square[ii][jj];
			}
			else
			{
				double r2 = (R - R0)*(R - R0) + Z*Z;
				psi_par = psi_0*(q_0 + 1) + 0.5*sqrt(4 * pow(psi_0*(q_0 + 1), 2.0) - 4 * psi_0*B0*(r2 - a*a));
			}
		}

		RHS_X = (1.0 / B_ss*(marker[p].v_par*B_s + times(b_par, alpha*mu / e*grad_B_par - deltaE_par)));
		RHS_V = (1.0 / m / B_ss*B_s*(e*deltaE_par / alpha - mu*grad_B_par));
		RHS_W = (-df0dt) / marker[p].g;


		if (SLOWING_DOWN)
		{
			//record old particle information
			v = sqrt(pow(marker[p].v_par, 2.0) + pow(marker[p].v_per, 2.0));
			lambda = marker[p].v_par / v;


			//slowing down particle process
			RHS_V += lambda*(-nu*(v + pow(v_c, 3.0) / pow(v, 2.0)));
			RHS_U = sqrt(1 - pow(lambda, 2.0))*(-nu*(v + pow(v_c, 3.0) / pow(v, 2.0)));

			//update particle information
			marker[p].v_per = Uold + timestep*RHS_U;
			Ua += coeff*dt*RHS_U;

			//mu should change
			mu = 0.5*species[0].mass*pow(marker[p].v_per, 2.0) / B_value_par;

		}

		RHS_X[2] = RHS_X[2] / marker[p].X[0];

		marker[p].X = Xold + timestep*RHS_X;
		marker[p].v_par = Vold + timestep*RHS_V;
		marker[p].w = Wold + timestep*RHS_W;
		Xa += coeff*dt*RHS_X;
		Va += coeff*dt*RHS_V;
		Wa += coeff*dt*RHS_W;
		t = told + timestep;


		//recycle particles
		if (RECYCLE_METHOD == 1)
		{
			if (NUMERICAL_EQUILIBRIUM)
			{
				if (abs(psi_par) < 1e-3)
				{
					OUT_OF_BOUDARY = true;
					break;
				}
			}
			else
			{
				if (sqrt(pow(marker[p].X[0] - R0, 2.0) + pow(marker[p].X[1], 2.0)) >= a)
				{
					OUT_OF_BOUDARY = true;
					break;
				}
			}
		}
		else if (RECYCLE_METHOD == 2)
		{
			if (NUMERICAL_EQUILIBRIUM)
			{
				if (abs(psi_par) < 1e-8)
				{
					marker[p].id = -1;
					break;
				}
			}
			else
			{
				if (sqrt(pow(marker[p].X[0] - R0, 2.0) + pow(marker[p].X[1], 2.0)) >= a)
				{
					marker[p].id = -1;
					break;
				}
			}
		}
	}


	if (OUT_OF_BOUDARY)
	{
		marker[p].v_par = Vold;
		marker[p].X = Xold;
		marker[p].w = Wold;
		marker[p].X[1] = -marker[p].X[1];
		if (SLOWING_DOWN)
			marker[p].v_per = Uold;
		//add by 2012-3-14
		if (marker[p].id == 0)
			LOSS = true;
	}
	else
	{
		marker[p].X = Xa;
		marker[p].v_par = Va;
		marker[p].w = Wa;
		if (SLOWING_DOWN)
			marker[p].v_per = Ua;
		//add by 2012-3-16
		if (TIMESTEP == 1)
			dZ_of_step1 = marker[p].X[1] - Xold[1];
		if (marker[p].X[1] * dZ_of_step1 < 0)
			Z_reverse_sign = true;
	}


	//full f
	if (FULL_F)
		marker[p].w = marker[p].f_over_g;


	i = (int)((marker[p].X[0] - R0 + a*a_b) / dR);
	j = (int)((marker[p].X[1] + a*a_b) / dZ);


	R = R0 - a*a_b + i*dR;
	Z = -a*a_b + j*dZ;



	weight_line_auxiliary[0][0] = (marker[p].X[0] - R) / dR;
	weight_line_auxiliary[1][0] = (marker[p].X[1] - Z) / dZ;

	weight_line_auxiliary[0][1] = 1.0 - weight_line_auxiliary[0][0];
	weight_line_auxiliary[1][1] = 1.0 - weight_line_auxiliary[1][0];


	for (int ii = 0; ii<2; ii++)
		for (int jj = 0; jj<2; jj++)
			weight_square[ii][jj] = weight_line_auxiliary[0][1 - ii] * weight_line_auxiliary[1][1 - jj];

	B_par = 0.0;
	curl_b_par = 0.0;
	grad_B_par = 0.0;
	B_value_par = 0.0;

	for (int ii = 0; ii<2; ii++)
		for (int jj = 0; jj<2; jj++)
		{
			B_value_par += B_value[i + ii][j + jj] * weight_square[ii][jj];
			for (int d = 0; d<2; d++)
			{
				B_par[d] += B[d][i + ii][j + jj] * weight_square[ii][jj];
				grad_B_par[d] += grad_B[d][i + ii][j + jj] * weight_square[ii][jj];
			}

			for (int d = 0; d<3; d++)
				curl_b_par[d] += curl_b[d][i + ii][j + jj] * weight_square[ii][jj];
		}

	if (NUMERICAL_EQUILIBRIUM)
	{
		double g_eq_par = 0.0;
		for (int ii = 0; ii<2; ii++)
			for (int jj = 0; jj<2; jj++)
				g_eq_par += g_eq[i + ii][j + jj] * weight_square[ii][jj];

		B_par[2] = g_eq_par / marker[p].X[0];
	}
	else
		B_par[2] = B0*R0 / marker[p].X[0];

	grad_B_par[2] = 0.0;

	for (int d = 0; d<3; d++)
		b_par[d] = B_par[d] / B_value_par;


	grad_deltaA_par_par = 0.0;
	deltaA_par_par = 0.0;
	deltaPhi_par = 0.0;

	for (int s = 0; s<NUM_MODE; s++)
	{
		grad_deltaPhi_X_par = 0.0;
		grad_deltaPhi_Y_par = 0.0;

		COS_NPHI_PLUS_OMEGA_T = cos(n[s] * marker[p].X[2] - omega_A[s] * t);
		SIN_NPHI_PLUS_OMEGA_T = sin(n[s] * marker[p].X[2] - omega_A[s] * t);

		for (int d = 0; d<3; d++)
			for (int I = i; I <= i + 1; I++)
				for (int J = j; J <= j + 1; J++)
				{

					if (d != 2)
					{
						grad_deltaPhi_X_local[d][I - i][J - j] = grad_deltaPhi_cos[s][d][I][J] * COS_NPHI_PLUS_OMEGA_T - grad_deltaPhi_sin[s][d][I][J] * SIN_NPHI_PLUS_OMEGA_T;
						grad_deltaPhi_Y_local[d][I - i][J - j] = grad_deltaPhi_sin[s][d][I][J] * COS_NPHI_PLUS_OMEGA_T + grad_deltaPhi_cos[s][d][I][J] * SIN_NPHI_PLUS_OMEGA_T;
					}
					else
					{
						grad_deltaPhi_X_local[d][I - i][J - j] = (-grad_deltaPhi_cos[s][d][I][J] * SIN_NPHI_PLUS_OMEGA_T - grad_deltaPhi_sin[s][d][I][J] * COS_NPHI_PLUS_OMEGA_T)*n[s] / marker[p].X[0];
						grad_deltaPhi_Y_local[d][I - i][J - j] = (-grad_deltaPhi_sin[s][d][I][J] * SIN_NPHI_PLUS_OMEGA_T + grad_deltaPhi_cos[s][d][I][J] * COS_NPHI_PLUS_OMEGA_T)*n[s] / marker[p].X[0];
					}
					if (d != 2)
						grad_deltaA_par_local[d][I - i][J - j] = X[s] * grad_deltaA_par_cos[s][d][I][J] * COS_NPHI_PLUS_OMEGA_T - X[s] * grad_deltaA_par_sin[s][d][I][J] * SIN_NPHI_PLUS_OMEGA_T + Y[s] * grad_deltaA_par_sin[s][d][I][J] * COS_NPHI_PLUS_OMEGA_T + Y[s] * grad_deltaA_par_cos[s][d][I][J] * SIN_NPHI_PLUS_OMEGA_T;
					else
						grad_deltaA_par_local[d][I - i][J - j] = (-X[s] * grad_deltaA_par_cos[s][d][I][J] * SIN_NPHI_PLUS_OMEGA_T - X[s] * grad_deltaA_par_sin[s][d][I][J] * COS_NPHI_PLUS_OMEGA_T - Y[s] * grad_deltaA_par_sin[s][d][I][J] * SIN_NPHI_PLUS_OMEGA_T + Y[s] * grad_deltaA_par_cos[s][d][I][J] * COS_NPHI_PLUS_OMEGA_T)*n[s] / marker[p].X[0];
				}

		for (int I = i; I <= i + 1; I++)
			for (int J = j; J <= j + 1; J++)
			{
				deltaA_par_local[I - i][J - j] = X[s] * deltaA_par_cos[s][I][J] * COS_NPHI_PLUS_OMEGA_T - X[s] * deltaA_par_sin[s][I][J] * SIN_NPHI_PLUS_OMEGA_T + Y[s] * deltaA_par_sin[s][I][J] * COS_NPHI_PLUS_OMEGA_T + Y[s] * deltaA_par_cos[s][I][J] * SIN_NPHI_PLUS_OMEGA_T;
				deltaPhi_local[I - i][J - j] = X[s] * deltaPhi_cos[s][I][J] * COS_NPHI_PLUS_OMEGA_T - X[s] * deltaPhi_sin[s][I][J] * SIN_NPHI_PLUS_OMEGA_T + Y[s] * deltaPhi_sin[s][I][J] * COS_NPHI_PLUS_OMEGA_T + Y[s] * deltaPhi_cos[s][I][J] * SIN_NPHI_PLUS_OMEGA_T;
			}


		for (int ii = 0; ii<2; ii++)
			for (int jj = 0; jj<2; jj++)
			{
				for (int d = 0; d<3; d++)
				{
					grad_deltaA_par_par[d] += grad_deltaA_par_local[d][ii][jj] * weight_square[ii][jj];
					grad_deltaPhi_X_par[d] += grad_deltaPhi_X_local[d][ii][jj] * weight_square[ii][jj];
					grad_deltaPhi_Y_par[d] += grad_deltaPhi_Y_local[d][ii][jj] * weight_square[ii][jj];

				}
				deltaA_par_par += deltaA_par_local[ii][jj] * weight_square[ii][jj];
				deltaPhi_par += deltaPhi_local[ii][jj] * weight_square[ii][jj];
			}

		deltaE_X_par[s] = (b_par*grad_deltaPhi_X_par)*b_par - grad_deltaPhi_X_par;
		deltaE_Y_par[s] = (b_par*grad_deltaPhi_Y_par)*b_par - grad_deltaPhi_Y_par;

	}
	deltaB_par = times(grad_deltaA_par_par, b_par) + deltaA_par_par*curl_b_par;

	if (SLOWING_DOWN)
		marker[p].mu = 0.5*species[0].mass*pow(marker[p].v_per, 2.0) / B_value_par;

	//update v_per
	marker[p].v_per = sqrt(2 * marker[p].mu*B_value_par / m);

	if (SCATTERING)
	{
		//if add collison effect,mu should be changed!!
		double v, lambda_old, lambda_new;
		double nu_d;
		const double c = 0.17*v_0;
		int pm;

		v = sqrt(pow(marker[p].v_par, 2.0) + pow(marker[p].v_per, 2.0));
		nu_d = nu*pow(v_c, 3.0) / 2.0 / (pow(c, 3.0) + pow(v, 3.0));

		lambda_old = marker[p].v_par / v;
		pm = (rand()*1.0 / RAND_MAX)>0.5 ? 1 : -1;
		lambda_new = lambda_old*(1 - 2 * nu_d*dt) + pm*sqrt((1 - lambda_old*lambda_old) * 2 * nu_d*dt);


		marker[p].v_par = v*lambda_new;
		marker[p].v_per = sqrt(v*v - pow(marker[p].v_par, 2.0));

		//mu should change
		marker[p].mu = 0.5*species[0].mass*pow(marker[p].v_per, 2.0) / B_value_par;
	}



	B_s = B_par + deltaB_par + alpha*m*marker[p].v_par / e*curl_b_par;
	B_ss = B_s*b_par;

	v_d = 1.0 / (e*B_ss)*(alpha*m*pow(marker[p].v_par, 2.0)*curl_b_par + marker[p].mu*times(b_par, grad_B_par));



	//update p_phi
	marker[p].P_phi = (m*marker[p].v_par*marker[p].X[0] * (B_par[2] / B_value_par) + e*psi_par / alpha);

	E = 0.5*m*pow(marker[p].v_par, 2.0) + marker[p].mu*B_value_par;

	//add by 2011-2.24 
	/*********************************************************************************************************
	*                                                                                                        *
	*                                   eliminate adiabatic term using                                       *
	*            0.retain adiabatic term : f_a = 0.0                                                         *
	*            1.candy's method        : f_a = e*R*deltaA_||*B_phi/B df_0dP_phi + e*deltaPhi*df_0dE        *
	*            2.fu's method           : f_a = -xi \cdot \grad f = e*R*deltaA_||*B_phi/B df_0dP_phi        *
	*                                                                                                        *
	*********************************************************************************************************/
	//adiabatic part of delta-f
	double deltaf_ad = 0.0;
	if (ADIABATIC != 0)
	{
		double psi_par;
		if (SHAFRANOV)
		{
			psi_par = 0;
			for (int ii = 0; ii<2; ii++)
				for (int jj = 0; jj<2; jj++)
					psi_par += psi[i + ii][j + jj] * weight_square[ii][jj];
		}
		else
		{
			R = marker[p].X[0];
			Z = marker[p].X[1];
			double r2 = (R - R0)*(R - R0) + Z*Z;
			psi_par = psi_0*(q_0 + 1) + 0.5*sqrt(4 * pow(psi_0*(q_0 + 1), 2.0) - 4 * psi_0*B0*(r2 - a*a));
		}

		v = sqrt(2 * E / species[0].mass);

		//p_phi =  (species[0].mass*marker[p].v_par*marker[p].X[0]*(B_par[2]/B_value_par)+species[0].charge*psi_par/alpha);
		p_phi = marker[p].P_phi;

		double bracket_psi;
		if (AVERAGE_PSI)
		{
			double mu = marker[p].mu;
			int sgn_v_parallel = marker[p].v_par>0 ? 1 : -1;
			if ((1 - mu*B0 / E) > 0)
				bracket_psi = p_phi / species[0].charge - species[0].mass / species[0].charge*sgn_v_parallel*v*R0*sqrt(1 - mu*B0 / E);
			else
				bracket_psi = p_phi / species[0].charge;
			double f0 = c_f*1.0 / (pow(v, 3) + v_c*v_c*v_c)*(1 + erf((v_0 - v) / (0.2*v_A)))*exp(-(bracket_psi / (c_1*Delta_psi)));

			df_0dP_phi = -1.0 / (e*c_1*(psi_a - psi_0))*f0;

			if ((1 - mu*B0 / E) > 0)
				df_0dE = -3 / m*sqrt(2 * E / m) / (pow(2 * E / m, 1.5) + pow(v_c, 3.0))*f0
				- 1 / (sqrt(2 * PI*E*m)*0.1*v_A)*exp(-(pow((v_0 - v) / (0.2*v_A), 2.0)))*c_f*1.0 / (pow(v, 3) + v_c*v_c*v_c)*exp(-(bracket_psi / (c_1*Delta_psi)))
				+ sgn_v_parallel*R0 / (species[0].charge*c_1*Delta_psi*sqrt(v*v - 2 * mu*B0 / species[0].mass))*f0;
			else
				df_0dE = -3 / m*sqrt(2 * E / m) / (pow(2 * E / m, 1.5) + pow(v_c, 3.0))*f0
				- 1 / (sqrt(2 * PI*E*m)*0.1*v_A)*exp(-(pow((v_0 - v) / (0.2*v_A), 2.0)))*c_f*1.0 / (pow(v, 3) + v_c*v_c*v_c)*exp(-(bracket_psi / (c_1*Delta_psi)));

		}
		else
		{
			df_0dP_phi = -1.0 / (e*c_1*Delta_psi)*c_f*1.0 / (pow(v, 3) + v_c*v_c*v_c)*(1 + erf((v_0 - v) / (0.2*v_A)))*exp(-(p_phi / (e*c_1*Delta_psi)));
			df_0dE = -3 / m*sqrt(2 * E / m) / (pow(2 * E / m, 1.5) + pow(v_c, 3.0))*c_f*1.0 / (pow(v, 3) + v_c*v_c*v_c)*(1 + erf((v_0 - v) / (0.2*v_A)))*exp(-(p_phi / (e*c_1*Delta_psi)))
				- 1 / (sqrt(2 * PI*E*m)*0.1*v_A)*exp(-(pow((v_0 - v) / (0.2*v_A), 2.0)))*c_f*1.0 / (pow(v, 3) + v_c*v_c*v_c)*exp(-(p_phi / (e*c_1*Delta_psi)));
		}
		if (ADIABATIC == 1)
			deltaf_ad = (e*marker[p].X[0] * deltaA_par_par*B_par[2] / B_value_par / alpha*df_0dP_phi + e*deltaPhi_par*df_0dE) / marker[p].g;
		else if (ADIABATIC == 2)
			deltaf_ad = (e*marker[p].X[0] * deltaA_par_par*B_par[2] / B_value_par / alpha*df_0dP_phi) / marker[p].g;

	}

	t = told;


	if (marker[p].id != -1)
	{
#pragma omp critical
		{
			for (int s = 0; s<NUM_MODE; s++)
			{
				myJ_dot_E_X[s] += e*(v_d + marker[p].v_par*deltaB_par / B_ss)*deltaE_X_par[s] * (marker[p].w - deltaf_ad);
				myJ_dot_E_Y[s] += e*(v_d + marker[p].v_par*deltaB_par / B_ss)*deltaE_Y_par[s] * (marker[p].w - deltaf_ad);
			}
		}
	}