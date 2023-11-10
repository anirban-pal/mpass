struct cBezier2 {
	//double PQB[48],PQB1[48];
  
  gsl_matrix * P = gsl_matrix_alloc (4, 3);
  gsl_matrix * Q = gsl_matrix_alloc (4, 3);

	gsl_matrix * PP = gsl_matrix_alloc (4, 4);
	gsl_matrix * QQ = gsl_matrix_alloc (4, 4);
  gsl_matrix * PQ = gsl_matrix_alloc (4, 4);
	
};

static int Integrand(const int *ndim, const double xx[],
  const int *ncomp, double ff[], void * params) {
  
  cBezier2 grparams = *(cBezier2 *) params; 
  
  double t1 = xx[0], t2 = xx[1];
  double t1c = 1.0-t1, t2c = 1.0-t2;

  double Xv[4] = {t1c*t1c*t1c, 3*t1c*t1c*t1, 3*t1c*t1*t1, t1*t1*t1};
  double Yv[4] = {t2c*t2c*t2c, 3*t2c*t2c*t2, 3*t2c*t2*t2, t2*t2*t2};

  double DXv[4] = {-3*t1c*t1c, 3*t1c*(t1c-2*t1), 3*t1*(2*t1c-t1), 3*t1*t1};
  double DYv[4] = {-3*t2c*t2c, 3*t2c*(t2c-2*t2), 3*t2*(2*t2c-t2), 3*t2*t2};
  
  double J1 = 0, J2 = 0, C = 0;

  for(int i=0;i<4;i++)
    for(int j=0;j<4;j++){
      J1 += DXv[i]*gsl_matrix_get (grparams.PP, i, j)*DXv[j];
      J2 += DYv[i]*gsl_matrix_get (grparams.QQ, i, j)*DYv[j];
    
      C += Xv[i]*gsl_matrix_get (grparams.PP, i, j)*Xv[j] 
          + Yv[i]*gsl_matrix_get (grparams.QQ, i, j)*Yv[j] 
          -2*Xv[i]*gsl_matrix_get (grparams.PQ, i, j)*Yv[j];
      
    }
  
  C = sqrt(C);
	ff[0] = V(C)*sqrt(J1)*sqrt(J2);	

  return 0;
}

static int I1234(const int *ndim, const double xx[],
  const int *ncomp, double F[], void * params) {
  
  cBezier2 grparams = *(cBezier2 *) params; 
  
  double t1 = xx[0], t2 = xx[1];
  double t1c = 1.0-t1, t2c = 1.0-t2;

  double Xv[4] = {t1c*t1c*t1c, 3*t1c*t1c*t1, 3*t1c*t1*t1, t1*t1*t1};
  double Yv[4] = {t2c*t2c*t2c, 3*t2c*t2c*t2, 3*t2c*t2*t2, t2*t2*t2};

  double DXv[4] = {-3*t1c*t1c, 3*t1c*(t1c-2*t1), 3*t1*(2*t1c-t1), 3*t1*t1};
  double DYv[4] = {-3*t2c*t2c, 3*t2c*(t2c-2*t2), 3*t2*(2*t2c-t2), 3*t2*t2};
  
  double J1 = 0, J2 = 0, C = 0;

  //double F[48] = {0};
  double cP[12] = {0},cQ[12] = {0},JP[12] = {0},JQ[12] = {0};

  for(int i=0;i<4;i++)
    for(int j=0;j<3;j++){
      //cP[i*3+j] = 0.0; cQ[i*3+j] = 0.0; JP[i*3+j] = 0.0; JQ[i*3+j] = 0.0;
      
      for(int k=0;k<4;k++) {
        cP[i*3+j] += Xv[i] * ( Xv[k] * gsl_matrix_get (grparams.P, k, j) - Yv[k] * gsl_matrix_get (grparams.Q, k, j) ) ;
        cQ[i*3+j] += Yv[i] * ( Yv[k] * gsl_matrix_get (grparams.Q, k, j) - Xv[k] * gsl_matrix_get (grparams.P, k, j) ) ;
        
        JP[i*3+j] += DXv[i] * DXv[k] * gsl_matrix_get (grparams.P, k, j) ;
        JQ[i*3+j] += DYv[i] * DYv[k] * gsl_matrix_get (grparams.Q, k, j) ;
        
      }
    }
  
  for(int i=0;i<4;i++)
    for(int j=0;j<4;j++){
      J1 += DXv[i]*gsl_matrix_get (grparams.PP, i, j)*DXv[j];
      J2 += DYv[i]*gsl_matrix_get (grparams.QQ, i, j)*DYv[j];
    
      C += Xv[i]*gsl_matrix_get (grparams.PP, i, j)*Xv[j] 
          + Yv[i]*gsl_matrix_get (grparams.QQ, i, j)*Yv[j] 
          -2*Xv[i]*gsl_matrix_get (grparams.PQ, i, j)*Yv[j];
      
    }
  
  C = sqrt(C);
  
  for(int i=0;i<12;i++)
  {
    F[i] = sqrt(J1*J2)*cP[i]*f(C);
    F[i+12] = sqrt(J2/J1)*JP[i]*V(C);
    F[i+24] = sqrt(J1*J2)*cQ[i]*f(C);
    F[i+36] = sqrt(J1/J2)*JQ[i]*V(C);
  }

  return 0;
}

double inter_energy(cBezier *b1, cBezier *b2, int force_flag)
{
  cBezier2 bpair;
  
  gsl_matrix_memcpy(bpair.P, b1->P);
  gsl_matrix_memcpy(bpair.Q, b2->P);
  
  gsl_matrix_memcpy(bpair.PP, b1->PP);
  gsl_matrix_memcpy(bpair.QQ, b2->PP);
  
  gsl_blas_dgemm (CblasNoTrans, CblasTrans, 1.0, bpair.P, bpair.Q, 0.0, bpair.PQ); 
  
  //gsl_matrix_print(bpair.PP,"PP");
  //gsl_matrix_print(bpair.QQ,"QQ");
  //gsl_matrix_print(bpair.PQ,"PQ");
  
  
  int comp, nregions, neval, nevalf, fail;
  double integral[1], error[1], prob[1];
  
  auto t1 = std::chrono::high_resolution_clock::now();
  
  Cuhre(2, 1, Integrand, &bpair, NVEC,
      EPSREL, EPSABS, VERBOSE | LAST,
      MINEVAL, MAXEVAL, KEY,
      STATEFILE, SPIN,
      &nregions, &neval, &fail, integral, error, prob);
  
  auto t2 = std::chrono::high_resolution_clock::now();
  
  double integralc[48], errorc[48], probc[48];    
  
  if (force_flag)
  Cuhre(2, 48, I1234, &bpair, NVEC,
      EPSREL, EPSABS, VERBOSE | LAST,
      MINEVAL, MAXEVAL, KEY,
      STATEFILE, SPIN,
      &nregions, &nevalf, &fail, integralc, errorc, probc);
  
  auto t3 = std::chrono::high_resolution_clock::now();
  
  b1->cohE += 0.5*integral[0];
  b2->cohE += 0.5*integral[0];
  
  loop(i,12) {
    b1->f[i] += (integralc[i]+integralc[i+12]);
    b2->f[i] += (integralc[i+24]+integralc[i+36]);
  }
  
  if(1) {
    auto dur1 = std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1);
    auto dur2 = std::chrono::duration_cast<std::chrono::nanoseconds>(t3 - t2);
    
    double E0 = -2.73525902052653568929940774978604e+00;
    double err = fabs((integral[0]-E0)/E0);
    printf("EPSREL: %le NEVAL %d Energy: %.32e Rel_error: %.32e Time (ns): %lld\n",EPSREL,neval,integral[0], err, dur1.count());
    
    loop(i,12) {
      //printf(" %.32e\n",integralc[i]+integralc[i+12]);
      //std::cout << integralc[i]+integralc[i+12] << std::endl;
      
    }
    
  }

  if(0) {
    double errorI = (integral[0]-b1->length0 * b2->length0)/(b1->length0 * b2->length0);
    printf("EPSREL=%le Length_product = %lf Integral = %lf Error = %le\n",EPSREL,b1->length0 * b2->length0, integral[0], errorI);
  } 

  return integral[0];
}

double inter_energy_discrete(cBezier *b1, cBezier *b2, int force_flag)
{
  cBezier2 bpair;
  
  gsl_matrix_memcpy(bpair.P, b1->P);
  gsl_matrix_memcpy(bpair.Q, b2->P);
  
  gsl_matrix_memcpy(bpair.PP, b1->PP);
  gsl_matrix_memcpy(bpair.QQ, b2->PP);
  
  gsl_blas_dgemm (CblasNoTrans, CblasTrans, 1.0, bpair.P, bpair.Q, 0.0, bpair.PQ); 
  
  auto t1 = std::chrono::high_resolution_clock::now();
  
  double F[48] = {0};
  double fac = (b1->length/b1->discrete.size())*(b2->length/b2->discrete.size());

  double c,U12 = 0;
  loop(i,b1->discrete.size())
    loop(j,b2->discrete.size()) {
      c = 0;
      loop2(k,1,4) {
        double d = b1->discrete[i].comp[k] - b2->discrete[j].comp[k]; c += d*d;
      }
      c = sqrt(c);
      
      U12 += V(c);    
      
      if(force_flag)
      {
        double t1 = b1->discrete[i].comp[0], t2 = b2->discrete[j].comp[0];
        double t1c = 1.0-t1, t2c = 1.0-t2;
      
        double Xv[4] = {t1c*t1c*t1c, 3*t1c*t1c*t1, 3*t1c*t1*t1, t1*t1*t1};
        double Yv[4] = {t2c*t2c*t2c, 3*t2c*t2c*t2, 3*t2c*t2*t2, t2*t2*t2};
      
        double DXv[4] = {-3*t1c*t1c, 3*t1c*(t1c-2*t1), 3*t1*(2*t1c-t1), 3*t1*t1};
        double DYv[4] = {-3*t2c*t2c, 3*t2c*(t2c-2*t2), 3*t2*(2*t2c-t2), 3*t2*t2};
        
        double J1 = 0, J2 = 0, C = 0;
      
        //double F[48] = {0};
        double cP[12] = {0},cQ[12] = {0},JP[12] = {0}, JQ[12] = {0};
      
        for(int ii=0;ii<4;ii++)
          for(int jj=0;jj<3;jj++){
            //cP[ii*3+jj] = 0.0; cQ[ii*3+jj] = 0.0; JP[ii*3+jj] = 0.0; JQ[ii*3+jj] = 0.0;
            
            for(int k=0;k<4;k++) {
              cP[ii*3+jj] += Xv[ii] * ( Xv[k] * gsl_matrix_get (b1->P, k, jj) - Yv[k] * gsl_matrix_get (b2->P, k, jj) ) ;
              cQ[ii*3+jj] += Yv[ii] * ( Yv[k] * gsl_matrix_get (b2->P, k, jj) - Xv[k] * gsl_matrix_get (b1->P, k, jj) ) ;
              
              JP[ii*3+jj] += DXv[ii] * DXv[k] * gsl_matrix_get (b1->P, k, jj) ;
              JQ[ii*3+jj] += DYv[ii] * DYv[k] * gsl_matrix_get (b2->P, k, jj) ;
              
            }
          }
        
        for(int ii=0;ii<4;ii++)
          for(int jj=0;jj<4;jj++){
            J1 += DXv[ii]*gsl_matrix_get (b1->PP, ii, jj)*DXv[jj];
            J2 += DYv[ii]*gsl_matrix_get (b2->PP, ii, jj)*DYv[jj];          
            
            //C += Xv[ii]*gsl_matrix_get (b1->PP, ii, jj)*Xv[jj] 
              //+ Yv[ii]*gsl_matrix_get (b2->PP, ii, jj)*Yv[jj] 
              //-2*Xv[ii]*gsl_matrix_get (bpair.PQ, ii, jj)*Yv[jj];
          }
        
        //C = sqrt(C);
        
        //printf("%lf %lf\n",c,C);
        for(int ii=0;ii<12;ii++)
        {
          F[ii] += fac*cP[ii]*f(c);
          F[ii+12] += (fac/J1)*JP[ii]*V(c);
          F[ii+24] += fac*cQ[ii]*f(c);
          F[ii+36] += (fac/J2)*JQ[ii]*V(c);
        }
      }
  }
  
  auto t2 = std::chrono::high_resolution_clock::now();

  loop(i,12) {
    b1->f[i] += (F[i]+F[i+12]);
    b2->f[i] += (F[i+24]+F[i+36]);
  }
  
  if(0)
  {
    loop(i,12) {
      printf("%e ",F[i]+F[i+12]);
      if((i+1)%4==0) printf("\n");
    }
    printf("\n");
    
    loop(i,12) {
      printf("%e ",F[i+24]+F[i+36]);
      if((i+1)%4==0) printf("\n");
    }
    printf("\n");
  }
  
  if(1) {
    auto dur1 = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1);
    //auto dur2 = std::chrono::duration_cast<std::chrono::microseconds>(t3 - t2);
  
    printf("\nN: %lu Energy: %.16e Time: %lld\n",b1->discrete.size(),U12*fac, dur1.count());
  }
  
  
  return U12*fac;
}


