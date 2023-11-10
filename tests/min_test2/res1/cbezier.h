class cBezier {
  
  public:
    double x[12],v[12],f[12],a[12],w[12];
    int CP[4], type, id, fiber;
    double energy, length, axialE, bendE, cohE, kinE, rho, mom[3];
    int leftB, rightB;
    double length0;
    std::vector<Vector4> discrete; 
    
    gsl_matrix * P = gsl_matrix_alloc (4, 3);
    gsl_matrix * PB = gsl_matrix_alloc (3, 4);
    gsl_matrix * P_B = gsl_matrix_alloc (4, 4);
		gsl_matrix * PP = gsl_matrix_alloc (4, 4);
    
    gsl_matrix * Q = gsl_matrix_alloc (4, 3);
    gsl_matrix * PQ_symm = gsl_matrix_alloc (4, 4);
    gsl_matrix * QQ = gsl_matrix_alloc (4, 4);
    gsl_matrix * M = gsl_matrix_alloc (4, 4);
    gsl_matrix * Mdot = gsl_matrix_alloc (4, 4);
    gsl_matrix * N = gsl_matrix_alloc (4, 4);
    
    double_t dM[4][4][4][3];
  
    void interpolateT(std::vector<Vector4> &, int);
    void interpolateL(std::vector<Vector4> &, int);
    void initialize_bezier();
    void update_bezier(int);
    void axial_ef();
    void bending_ef();
    void penalty_ef();
    void momentum();
    void kin_energy();
    
  private:
    void length_bezier();
    void set_length0();
    
};

void cBezier::interpolateT(std::vector<Vector4> &points, int n) {
  
  for(double t=1.0/n; t<(1.0+1.0/n); t+=1.0/n)
  {
    double tc = 1.0-t;
    double X[4] = {tc*tc*tc, 3*tc*tc*t, 3*tc*t*t, t*t*t};
    
    double rx = X[0]*x[0] + X[1]*x[3] + X[2]*x[6] + X[3]*x[9];
    double ry = X[0]*x[1] + X[1]*x[4] + X[2]*x[7] + X[3]*x[10];
    double rz = X[0]*x[2] + X[1]*x[5] + X[2]*x[8] + X[3]*x[11];
    
    points.push_back(Vector4(fiber,rx,ry,rz));
  }
}

double ds(double t, void *p)
{
	cBezier params = *(cBezier *) p;
	
  double tc = 1.0-t;
  double DX[4] = {-3*tc*tc, 3*tc*(tc-2*t), -3*t*(t-2*tc), 3*t*t};
  
	double J = 0;
	
  loop(i,4)
		loop(j,4)
			J += DX[i]*gsl_matrix_get(params.PP, i, j)*DX[j];
			
	return sqrt(fabs(J));	
}

void cBezier::interpolateL(std::vector<Vector4> &points, int n) {
  
  double result, error; size_t nev;
  
  gsl_function B1;  B1.function = &ds;  B1.params = this;
  gsl_integration_cquad_workspace * w1 = gsl_integration_cquad_workspace_alloc(100);

  //std::map<double, double> TLmap;
  size_t N = 10*n;
  
  //std::vector<double> lgrid, tgrid;
  double *lgrid = new double[N+1];
  double *tgrid = new double[N+1];
  
  int c=0;
  
  //gsl_matrix_print(this->P,"P");
  //gsl_matrix_print(this->PP,"PP");
  
  //for(double t=0.0; t<(1.0+0.1/n); t+=0.1/n)
  for(int i=0; i<=N; i++)
  {
    double t = i*1.0/N;
    gsl_integration_cquad (&B1, 0.0, t, ABS_ERR, REL_ERR2,w1, &result, &error, &nev);
    //TLmap.insert(std::make_pair(t, result));
    
    //tgrid.push_back(t); lgrid.push_back(result); c++;
    tgrid[c] = t; lgrid[c] = result;  c++;
    //printf("%d %lf %lf\n",this->id,t,result);
    
  }
  
  gsl_integration_cquad_workspace_free (w1);
  
  gsl_interp_accel *acc = gsl_interp_accel_alloc();
  gsl_spline *spline_steffen = gsl_spline_alloc(gsl_interp_steffen, N+1);
  gsl_spline_init(spline_steffen, lgrid, tgrid, N+1);

  int Ni = n; discrete.clear();
  for (int i = 1; i <= Ni; ++i)
  //for (float i = 0.5; i < Ni; ++i)
  {
    double li = (1 - i*1.0 / Ni) * lgrid[0] + (i*1.0 / Ni) * lgrid[N-1];
    double t = gsl_spline_eval(spline_steffen, li, acc);

    //printf("%g : %g\n", li, t);
    
    double tc = 1.0-t;
    double X[4] = {tc*tc*tc, 3*tc*tc*t, 3*tc*t*t, t*t*t};
    
    double rx = X[0]*x[0] + X[1]*x[3] + X[2]*x[6] + X[3]*x[9];
    double ry = X[0]*x[1] + X[1]*x[4] + X[2]*x[7] + X[3]*x[10];
    double rz = X[0]*x[2] + X[1]*x[5] + X[2]*x[8] + X[3]*x[11];
    
    points.push_back(Vector4(fiber,rx,ry,rz));
    discrete.push_back(Vector4(t,rx,ry,rz));
  }

  gsl_spline_free(spline_steffen);
  gsl_interp_accel_free(acc);

  delete[] lgrid;
  delete[] tgrid;
  
}

void cBezier::initialize_bezier() 
{
  double Bvec[] = {-1, 3, -3, 1, 3, -6, 3, 0, -3, 3, 0, 0, 1, 0, 0, 0};
  
  const gsl_matrix_const_view P1 = gsl_matrix_const_view_array( x, 4, 3 );
	const gsl_matrix_const_view B1 = gsl_matrix_const_view_array( Bvec, 4, 4 );
  
  gsl_matrix_memcpy(P, &P1.matrix);
  gsl_blas_dgemm (CblasTrans, CblasNoTrans, 1.0, &P1.matrix, &B1.matrix, 0.0, PB);
  gsl_blas_dgemm (CblasTrans, CblasNoTrans, 1.0, PB, PB, 0.0, P_B);
  gsl_blas_dgemm (CblasNoTrans, CblasTrans, 1.0, &P1.matrix, &P1.matrix, 0.0, PP);
  
  //gsl_matrix_print(PP,"PP:");
  if(fabs(length0)<1e-6) {
     length_bezier(); set_length0();
  }
  //printf("Initial Length = %lf\n",length);
  
  loop(i,12) {v[i] = 0.0, f[i] = 0.0;}
  
  update_bezier(0);
}

static int XXJ(const int *ndim, const double xx[],
  const int *ncomp, double ff[], void *p) {
  
  cBezier params = *(cBezier *) p;

  double t = xx[0], tc = 1.0 - xx[0], Jsq = 0;
  double X[4] = {tc*tc*tc, 3*tc*tc*t, 3*tc*t*t, t*t*t};
  double DX[4] = {-3*tc*tc, 3*tc*(tc-2*t), -3*t*(t-2*tc), 3*t*t};
  
  loop(i,4) loop(j,4) Jsq += DX[i]*gsl_matrix_get(params.PP, i, j)*DX[j];
  double J = sqrt(fabs(Jsq));
  
  loop(i,4) loop(j,4) {
    ff[i*4+j] = X[i]*X[j]*J;
    
    //printf("XXJ ffij: %lf - %lf %lf %lf\n",ff[i*4+j],X[i],X[j],J);
  }
  
  return 0;
}

static int XXf(const int *ndim, const double xx[],
  const int *ncomp, double ff[], void *p) {
  
  cBezier params = *(cBezier *) p;

  double t = xx[0], tc = 1.0 - xx[0], Jsq = 0, f = 0;
  double X[4] = {tc*tc*tc, 3*tc*tc*t, 3*tc*t*t, t*t*t};
  double DX[4] = {-3*tc*tc, 3*tc*(tc-2*t), -3*t*(t-2*tc), 3*t*t};
  
  loop(i,4) loop(j,4) Jsq += DX[i]*gsl_matrix_get(params.PP, i, j)*DX[j];
  double J = sqrt(fabs(Jsq));
  
  loop(i,4) loop(j,4) f += DX[i]*gsl_matrix_get(params.PQ_symm, i, j)*DX[j];
  //f = f/2;
  
  loop(i,4) loop(j,4) {
    ff[i*4+j] = X[i]*X[j]*f/J;
  }
  
  return 0;
}

static int YYf(const int *ndim, const double xx[],
  const int *ncomp, double ff[], void *p) {
  
  cBezier params = *(cBezier *) p;

  double t = xx[0], tc = 1.0 - xx[0], Jsq = 0, f = 0;
  double X[4] = {tc*tc*tc, 3*tc*tc*t, 3*tc*t*t, t*t*t};
  double DX[4] = {-3*tc*tc, 3*tc*(tc-2*t), -3*t*(t-2*tc), 3*t*t};
  
  loop(i,4) loop(j,4) Jsq += DX[i]*gsl_matrix_get(params.PP, i, j)*DX[j];
  double J = sqrt(fabs(Jsq));
  
  loop(i,4) loop(j,4) f += X[i]*gsl_matrix_get(params.QQ, i, j)*X[j];
  f = f/2;
  
  loop(i,4) loop(j,4) {
    ff[i*4+j] = (params.rho/2.0)*DX[i]*DX[j]*f/J - (params.kinE/params.length)*DX[i]*DX[j]/J;
  }
  
  return 0;
}

static int XXYYQ(const int *ndim, const double xx[],
  const int *ncomp, double ff[], void *p) {
  
  cBezier params = *(cBezier *) p;

  double t = xx[0], tc = 1.0 - xx[0], Jsq = 0, f = 0;
  double X[4] = {tc*tc*tc, 3*tc*tc*t, 3*tc*t*t, t*t*t};
  double DX[4] = {-3*tc*tc, 3*tc*(tc-2*t), -3*t*(t-2*tc), 3*t*t};
  
  loop(i,4) loop(j,4) Jsq += DX[i]*gsl_matrix_get(params.PP, i, j)*DX[j];
  double J = sqrt(fabs(Jsq));
  
  loop(i,4) loop(j,4) loop(m,4) loop(n,3) 
  {
    ff[i*4*4*3 + j*4*3 + m*3 + n] = 0;
    loop(k,4) ff[i*4*4*3 + j*4*3 + m*3 + n] += X[i]*X[j]*DX[m]*DX[k]*gsl_matrix_get(params.P, k, n)/J;
  }
  
  return 0;
}

static int YYQ(const int *ndim, const double xx[],
  const int *ncomp, double ff[], void *p) {
  
  cBezier params = *(cBezier *) p;

  double t = xx[0], tc = 1.0 - xx[0], Jsq = 0, f = 0;
  double X[4] = {tc*tc*tc, 3*tc*tc*t, 3*tc*t*t, t*t*t};
  double DX[4] = {-3*tc*tc, 3*tc*(tc-2*t), -3*t*(t-2*tc), 3*t*t};
  
  loop(i,4) loop(j,4) Jsq += DX[i]*gsl_matrix_get(params.PP, i, j)*DX[j];
  double J = sqrt(fabs(Jsq));
  
  loop(m,4) loop(n,3) 
  {
    ff[m*3 + n] = 0;
    loop(k,4) ff[m*3 + n] += DX[m]*DX[k]*gsl_matrix_get(params.P, k, n)/J;
  }
  
  return 0;
}

void cBezier::update_bezier(int dm_flag) 
{
  double Bvec[] = {-1, 3, -3, 1, 3, -6, 3, 0, -3, 3, 0, 0, 1, 0, 0, 0};
  
  const gsl_matrix_const_view P1 = gsl_matrix_const_view_array( x, 4, 3 );
	const gsl_matrix_const_view B1 = gsl_matrix_const_view_array( Bvec, 4, 4 );
  
  gsl_matrix_memcpy(P, &P1.matrix);
  gsl_blas_dgemm (CblasTrans, CblasNoTrans, 1.0, &P1.matrix, &B1.matrix, 0.0, PB);
  gsl_blas_dgemm (CblasTrans, CblasNoTrans, 1.0, PB, PB, 0.0, P_B);
  gsl_blas_dgemm (CblasNoTrans, CblasTrans, 1.0, &P1.matrix, &P1.matrix, 0.0, PP);
  
  
  length_bezier();
  rho = rho0[type-1]*(length0/length);
  kin_energy();
  
  //gsl_matrix_print(PP,"bezier PP:");
  
  //update M and Mdot
  int comp, nregions, neval, fail;
  double integral0[16], error0[16], prob0[16];
  
  Cuhre(2, 16, XXJ, this, NVEC,
    EPSREL, EPSABS, VERBOSE | LAST,
    MINEVAL, MAXEVAL, KEY,
    STATEFILE, SPIN,
    &nregions, &neval, &fail, integral0, error0, prob0);
  
  loop(i,4) loop(j,4) {
    //printf("Mij: %lf\n",integral0[i*4+j]);
    gsl_matrix_set(M,i,j,rho*integral0[i*4+j]);
  }
  
  //gsl_matrix_print(M,"bezier M:");
  
  if(dm_flag==0) {
  
    const gsl_matrix_const_view Q1 = gsl_matrix_const_view_array( v, 4, 3 );
    gsl_matrix_memcpy(Q, &Q1.matrix);
    
    //gsl_blas_dgemm (CblasNoTrans, CblasTrans, 1.0, Q, P, 0.0, PQ_symm);
    gsl_blas_dgemm (CblasNoTrans, CblasTrans, 1.0, P, Q, 0.0, PQ_symm);
    
    double integral[16], error[16], prob[16];
    
    Cuhre(2, 16, XXf, this, NVEC,
      EPSREL, EPSABS, VERBOSE | LAST,
      MINEVAL, MAXEVAL, KEY,
      STATEFILE, SPIN,
      &nregions, &neval, &fail, integral, error, prob);
      
    double fsum = 0.0;
    loop(i,4) loop(j,4) fsum += integral[i*4+j];
    
    loop(i,4) loop(j,4) {
      double val = rho*integral[i*4+j] - (fsum/length)*gsl_matrix_get(M,i,j);
      gsl_matrix_set(Mdot,i,j,val);  
    }
    
    gsl_blas_dgemm (CblasNoTrans, CblasTrans, 1.0, Q, Q, 0.0, QQ);
    
    Cuhre(2, 16, YYf, this, NVEC,
      EPSREL, EPSABS, VERBOSE | LAST,
      MINEVAL, MAXEVAL, KEY,
      STATEFILE, SPIN,
      &nregions, &neval, &fail, integral, error, prob);
      
    loop(i,4) loop(j,4) gsl_matrix_set(N,i,j,integral[i*4+j]);  
  }
  //gsl_matrix_print(PP,"PP:");
  //length_bezier(); set_length0();
  //printf("Length = %lf\n",length);
  
  loop(i,12) {v[i] = 0.0, f[i] = 0.0;}
  
  if(dm_flag==1) {
    //update dMdQ
    //printf("\nComputing dMdQ\n");
    int comp, nregions, neval, fail;
    double integral[192], error[192], prob[192];
    
    Cuhre(2, 192, XXYYQ, this, NVEC,
      EPSREL, EPSABS, VERBOSE | LAST,
      MINEVAL, MAXEVAL, KEY,
      STATEFILE, SPIN,
      &nregions, &neval, &fail, integral, error, prob);
    
    double integral2[12], error2[12], prob2[12];
    
    Cuhre(2, 12, YYQ, this, NVEC,
      EPSREL, EPSABS, VERBOSE | LAST,
      MINEVAL, MAXEVAL, KEY,
      STATEFILE, SPIN,
      &nregions, &neval, &fail, integral2, error2, prob2);
    
    loop(i,4) loop(j,4) loop(m,4) loop(n,3) 
      dM[i][j][m][n] = rho*integral[i*4*4*3 + j*4*3 + m*3 + n] - gsl_matrix_get(M,i,j)*integral2[m*3 + n]/length;
    
    //loop(i,4) loop(j,4) gsl_matrix_set(M,i,j,rho*integral[i*4+j]);
  
  }
}

void cBezier::length_bezier()
{
  double result, error; size_t nev;
  gsl_function B1;
  B1.function = &ds;
  B1.params = this;
  gsl_integration_cquad_workspace * w1 = gsl_integration_cquad_workspace_alloc(100);
  
  gsl_integration_cquad (&B1, 0.0, 1.0, ABS_ERR, REL_ERR2,w1, &result, &error, &nev);
  length = result;
  gsl_integration_cquad_workspace_free (w1);
}

void cBezier::set_length0() 
{
  length0 = length;
}

static int dJ(const int *ndim, const double xx[],
  const int *ncomp, double ff[], void *p) {
  
  cBezier params = *(cBezier *) p;

  double t = xx[0], tc = 1.0 - xx[0], Jsq = 0;
  double DX[4] = {-3*tc*tc, 3*tc*(tc-2*t), -3*t*(t-2*tc), 3*t*t};
  
  loop(i,4) loop(j,4) Jsq += DX[i]*gsl_matrix_get(params.PP, i, j)*DX[j];
  
  loop(i,4) loop(j,3) {
    double sum = 0.0;
    loop(k,4) sum += DX[i]*DX[k]*gsl_matrix_get(params.P, k, j);
    ff[i*3+j] = sum/sqrt(fabs(Jsq));
  }
  
  return 0;
}

void cBezier::axial_ef()
{
  length_bezier();
  axialE = 0.5*EA[type-1]*(length-length0)*(length-length0);
  
  //printf("\n Axial Energy: %lf length %lf",axialE,length);
  
  int comp, nregions, neval, fail;
  double integral[12], error[12], prob[12];
  
  Cuhre(2, 12, dJ, this, NVEC,
    EPSREL, EPSABS, VERBOSE | LAST,
    MINEVAL, MAXEVAL, KEY,
    STATEFILE, SPIN,
    &nregions, &neval, &fail, integral, error, prob);
    
  loop(i,12) {
    f[i] += -0.5*EA[type-1]*(length-length0)*integral[i];
    //f[i] += integral[i];
  }
}

void cBezier::penalty_ef()
{
  double x01[3],x32[3],x12[3],l01,l12,l32;
  
  loop(i,3) {
    x01[i] = gsl_matrix_get(P, 1, i) - gsl_matrix_get(P, 0, i);
    x32[i] = gsl_matrix_get(P, 2, i) - gsl_matrix_get(P, 3, i);
    x12[i] = gsl_matrix_get(P, 2, i) - gsl_matrix_get(P, 1, i);	    

    l01 += x01[i]*x01[i]; l12 += x12[i]*x12[i]; l32 += x32[i]*x32[i];
  } 
  
  l01 = sqrt(l01); l12 = sqrt(l12); l32 = sqrt(l32);
 
  if(0) {
    if(l01 > 0.5*length0) {
      loop(i,3) {
        f[3+i] += -2*PEN_STIFF*x01[i];
        f[i] += 2*PEN_STIFF*x01[i];
      }
    }
  
    if(l01 < 0.1*length0) {
      loop(i,3) {
        f[3+i] += 2*PEN_STIFF*x01[i];
        f[i] += -2*PEN_STIFF*x01[i];
      }
    }
    
    if(l32 > 0.5*length0) {
      loop(i,3) {
        f[6+i] += -2*PEN_STIFF*x32[i];
        f[9+i] += 2*PEN_STIFF*x32[i];
      }
    }
  
    if(l32 < 0.1*length0) {
      loop(i,3) {
        f[6+i] += 2*PEN_STIFF*x32[i];
        f[9+i] += -2*PEN_STIFF*x32[i];
      }
    }
  }
  
  //printf("\n Axial Energy: %lf length %lf",axialE,length);
}

double kappa_squared_ds(double t, void *p)
{
  cBezier params = *(cBezier *) p;
	
  double tc = 1.0-t;
  
  double X[4] = {tc*tc*tc, 3*tc*tc*t, 3*tc*t*t, t*t*t};
  double Y[4] = {-3*tc*tc, 3*tc*(tc-2*t), -3*t*(t-2*tc), 3*t*t};
  double Z[4] = {6*tc, 3*(2-6*tc), 3*(2-6*t), 6*t};
    
  double dr2=0, ddr2=0, drddr=0;
  loop(i,4) loop(j,4) {
    dr2 += Y[i]*gsl_matrix_get(params.PP, i, j)*Y[j];
    ddr2 += Z[i]*gsl_matrix_get(params.PP, i, j)*Z[j];
    drddr += Y[i]*gsl_matrix_get(params.PP, i, j)*Z[j];
  }
  
  double drxddr2 = dr2*ddr2-drddr*drddr;
  double kappa2 = drxddr2/(dr2*dr2*dr2);
  
  if(fabs(dr2)<=DR2_ERR) kappa2=0.0; //This appears reasonable and correct.
  
  //gsl_matrix_print(params.P,"params.P");
  //gsl_matrix_print(params.PP,"params.PP");
  //printf("kappa2ds %lf %.16f %lf %lf %.16f %.16f\n",t,dr2,ddr2,drddr,drxddr2,kappa2);
  
	return kappa2*sqrt(fabs(dr2));	
}

static int dB(const int *ndim, const double xx[],
  const int *ncomp, double ff[], void *p) {
  
  cBezier params = *(cBezier *) p;

  double t = xx[0], tc = 1.0 - xx[0], Jsq = 0;
  
  double X[4] = {tc*tc*tc, 3*tc*tc*t, 3*tc*t*t, t*t*t};
  double Y[4] = {-3*tc*tc, 3*tc*(tc-2*t), -3*t*(t-2*tc), 3*t*t};
  double Z[4] = {6*tc, 3*(2-6*tc), 3*(2-6*t), 6*t};
    
  double dr2=0, ddr2=0, drddr=0;
  loop(i,4) loop(j,4) {
    dr2 += Y[i]*gsl_matrix_get(params.PP, i, j)*Y[j]; // |r'|^2
    ddr2 += Z[i]*gsl_matrix_get(params.PP, i, j)*Z[j]; // |r''|^2
    drddr += Y[i]*gsl_matrix_get(params.PP, i, j)*Z[j]; // |r'.r''|
  }
  
  double drxddr2 = dr2*ddr2-drddr*drddr; 
  double kappa2 = drxddr2/(dr2*dr2*dr2);
  
  if(fabs(dr2)<=DR2_ERR ) kappa2=0.0; //This appears reasonable and correct.
  
  loop(i,4) loop(j,3) {
    double sum1 = 0.0, sum2 = 0.0;
    loop(k,4) { 
      sum1 += 2*dr2*dr2*dr2*(  dr2*Z[i]*Z[k] + ddr2*Y[i]*Y[k] - drddr*(Y[i]*Z[k]+Z[i]*Y[k]) )*gsl_matrix_get(params.P, k, j);
      sum1 -= 6*drxddr2*dr2*dr2*Y[i]*Y[k]*gsl_matrix_get(params.P, k, j);
      sum2 += Y[i]*Y[k]*gsl_matrix_get(params.P, k, j);
    }
    
    //double ffij0 = (sum1*dr2/pow(dr2,6) + sum2*kappa2)/sqrt(dr2);
    
    double term1 = sum1*sqrt(fabs(dr2))/pow(dr2,6);
    double term2 = sum2*kappa2/sqrt(fabs(dr2));
    
    if(fabs(sum1)<=KAPPA_ERR) term1 = 0.0; //Check again later
    if(fabs(sum2)<=KAPPA_ERR) term2 = 0.0; //Check again later
    
    double ffij = term1+term2;
    //printf("kappa2ds %lf %lf %lf %lf %.12f %lf\n",t,dr2,ddr2,drddr,drxddr2,kappa2); 
    //printf("ffij %lf %lf %lf %lf %lf %lf %lf\n",dr2,kappa2,sum1,sum2,term1,term2,ffij);
    
    ff[i*3+j] = ffij;
  }
  
  return 0;
}

void cBezier::bending_ef()
{
  double resultc, errorc; size_t nev;
  gsl_function B1;
  B1.function = &kappa_squared_ds;
  B1.params = this;
  gsl_integration_cquad_workspace * w1 = gsl_integration_cquad_workspace_alloc(100);
  
  gsl_integration_cquad (&B1, 0.0, 1.0, ABS_ERR, REL_ERR2,w1, &resultc, &errorc, &nev);
  
  bendE = 0.5*EI[type-1]*resultc;

  if(0) printf("nbE %16.6f\n",resultc);

  if(0) {
    
    Vector3 c03(x[9]-x[0],x[10]-x[1],x[11]-x[2]);
    Vector3 c01(x[3]-x[0],x[4]-x[1],x[5]-x[2]);
    Vector3 c23(x[9]-x[6],x[10]-x[7],x[11]-x[8]);
    
    double c03_1 = cross_mag(c01,c03)/(vec_mag(c01)*vec_mag(c03));
    double c03_2 = cross_mag(c23,c03)/(vec_mag(c23)*vec_mag(c03));
    
    printf("nbE %16.6f \t Ltest %16.6f %16.6f\n",resultc,c03_1,c03_2);
    
  }

  int comp, nregions, neval, fail;
  double integral[12], error[12], prob[12];
  
  Cuhre(2, 12, dB, this, NVEC,
    EPSREL, EPSABS, VERBOSE | LAST,
    MINEVAL, MAXEVAL, KEY,
    STATEFILE, SPIN,
    &nregions, &neval, &fail, integral, error, prob);
    
  loop(i,12) {
    f[i] += -0.5*EI[type-1]*integral[i];
  }

  if(0) loop(i,4) {
    printf("\n%16.6f %16.6f %16.6f",integral[i*3],integral[i*3+1],integral[i*3+2]);
  }
  
  //printf("%lf : %lf %lf %lf\n",resultc,integral[0],integral[1],integral[2]);
  
  gsl_integration_cquad_workspace_free (w1);
}

static int XJ(const int *ndim, const double xx[],
  const int *ncomp, double ff[], void *p) {
  
  cBezier params = *(cBezier *) p;

  double t = xx[0], tc = 1.0 - xx[0], Jsq = 0;
  double X[4] = {tc*tc*tc, 3*tc*tc*t, 3*tc*t*t, t*t*t};
  double DX[4] = {-3*tc*tc, 3*tc*(tc-2*t), -3*t*(t-2*tc), 3*t*t};
  
  loop(i,4) loop(j,4) Jsq += DX[i]*gsl_matrix_get(params.PP, i, j)*DX[j];
  double J = sqrt(fabs(Jsq));
  
  loop(i,4) {
    ff[i] = X[i]*J;
  }
  
  return 0;
}

void cBezier::momentum()
{
  int comp, nregions, neval, fail;
  double integral[4], error[4], prob[4];
  
  Cuhre(2, 4, XJ, this, NVEC,
    EPSREL, EPSABS, VERBOSE | LAST,
    MINEVAL, MAXEVAL, KEY,
    STATEFILE, SPIN,
    &nregions, &neval, &fail, integral, error, prob);
  
  const gsl_matrix_const_view Q1 = gsl_matrix_const_view_array( v, 4, 3 );
  gsl_matrix_memcpy(Q, &Q1.matrix);  
  
  //gsl_matrix_print(Q,"Q");
  
  loop(i,3) {
    mom[i] = 0;
    loop(j,4) {
      mom[i] += rho*gsl_matrix_get(Q,j,i)*integral[j];
      //printf("mom[%d] += %.12f = %.12f\n",i,rho*gsl_matrix_get(Q,j,i)*integral[j], mom[i]);
    }
    //printf("mom[%d]: %.12f\n",i,mom[i]);
  }
  
  //printf("rho %lf integral %lf %lf %lf %lf Bezier mom: %lf %lf %lf\n",rho,integral[0],integral[1],integral[2],integral[3],mom[0],mom[1],mom[2]);
}

double vel_squared_ds(double t, void *p)
{
	cBezier params = *(cBezier *) p;
	
  double tc = 1.0-t, Jsq = 0;
    
  double X[4] = {tc*tc*tc, 3*tc*tc*t, 3*tc*t*t, t*t*t};
  double DX[4] = {-3*tc*tc, 3*tc*(tc-2*t), -3*t*(t-2*tc), 3*t*t};
  
  double vv=0;
  loop(i,4) loop(j,4) {
    Jsq += DX[i]*gsl_matrix_get(params.PP, i, j)*DX[j];
    vv += X[i]*gsl_matrix_get(params.QQ, i, j)*X[j];
  }
  
	return vv*sqrt(fabs(Jsq));	
}

void cBezier::kin_energy()
{
  double resultc, errorc; size_t nev;
  
  gsl_function B1;
  B1.function = &vel_squared_ds;
  B1.params = this;
  
  gsl_integration_cquad_workspace * w1 = gsl_integration_cquad_workspace_alloc(100);
  gsl_integration_cquad (&B1, 0.0, 1.0, ABS_ERR, REL_ERR2,w1, &resultc, &errorc, &nev);
  
  kinE = 0.5*rho0[type-1]*(length0/length)*resultc;
  
  gsl_integration_cquad_workspace_free (w1);
}
