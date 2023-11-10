void print_QPXY(double *Q, double *P, double *X, double *Y, int np)
{
  loop(i,np) {
    printf("\n%10.2f %10.2f %10.2f\t%10.2f %10.2f %10.2f\t%10.2f %10.2f %10.2f\t%10.2f %10.2f %10.2f",
	Q[i*3],Q[i*3+1],Q[i*3+2],P[i*3],P[i*3+1],P[i*3+2],X[i*3],X[i*3+1],X[i*3+2],Y[i*3],Y[i*3+1],Y[i*3+2]);

  }	
}

void dH(const double q[], const double p[], void *fn_data, int cflag)
{
  fibnetwork *FN = (fibnetwork *) fn_data;

  int np = FN->np;	int nb = FN->nb;  int nf = FN->nf;
  
  //printf("\nq[]:\n");
  //loop(i,np) {
    //loop(j,3) {
      //printf(" %lf",q[i*3+j]);
    //}
    //printf("\n");
  //}
  
  //printf("\np[]:\n");
  //loop(i,np) {
    //loop(j,3) {
      //printf(" %lf",p[i*3+j]);
    //}
    //printf("\n");
  //}
    
  
  //printf("\n");
  //loop(i,np) printf("funcSTEP1: %d %lf %lf %lf %lf %lf %lf\n",i+1,FN->CPx[i].comp[0],FN->CPx[i].comp[1],FN->CPx[i].comp[2],FN->CPv[i].comp[0],FN->CPv[i].comp[1],FN->CPv[i].comp[2]);
  
  loop(i,np) loop(j,3) if(!isnan(q[i*3+j])) FN->CPx[i].comp[j] = q[i*3+j];
  loop(i,np) loop(j,3) if(!isnan(p[i*3+j])) FN->CPw[i].comp[j] = p[i*3+j]; //define momenta from momenta vector, define velocity later
  
  //apply fixed constraints
  if (1) 
  for(std::map<int,double>::iterator it = FN->fixed_dof.begin(); it != FN->fixed_dof.end(); it++) {
    int index = it->first-1;
    int i0 = index/3;
    int j0 = index%3;
    FN->CPx[i0].comp[j0] = it->second;
    FN->CPv[i0].comp[j0] = 0.0;
    //printf("FC APPLY: %d %d %d %lf %lf\n",i0,j0,index,it->second,FN->CPx[i0].comp[j0]);
  }
  
  //apply C1/C2
  if(cflag)
    loop(i,FN->C2_cps.size()) {
      int a = FN->C2_cps[i].comp[0];
      int b = FN->C2_cps[i].comp[1];
      int c = FN->C2_cps[i].comp[2];
      int d = FN->C2_cps[i].comp[3];
      int e = FN->C2_cps[i].comp[4];
      
      Vector3 db(FN->CPx[d-1].comp[0]-FN->CPx[b-1].comp[0], FN->CPx[d-1].comp[1]-FN->CPx[b-1].comp[1] , FN->CPx[d-1].comp[2]-FN->CPx[b-1].comp[2] );
      Vector3 eb(FN->CPx[e-1].comp[0]-FN->CPx[b-1].comp[0], FN->CPx[e-1].comp[1]-FN->CPx[b-1].comp[1] , FN->CPx[e-1].comp[2]-FN->CPx[b-1].comp[2] );
      Vector3 ad(FN->CPx[a-1].comp[0]-FN->CPx[d-1].comp[0], FN->CPx[a-1].comp[1]-FN->CPx[d-1].comp[1] , FN->CPx[a-1].comp[2]-FN->CPx[d-1].comp[2] );
      
      double m = sqrt(cross_mag(db,eb)/cross_mag(db,ad));
      double x = 1.0/(1.0+m);
      
      //EDIT to avoid singularities from control point overlap
      if(x<0.05) x = 0.05; if(x>0.95) x = 0.95;
      
      //x = 0.5;
      
      //printf("%d %d %d %d %d %d %lf %lf\n",i,a,b,c,d,e,m,x);
      loop(j,3) FN->CPx[c-1].comp[j] = x*FN->CPx[d-1].comp[j] + (1.0-x)*FN->CPx[b-1].comp[j];
    }
  
  //loop(i,np) printf("funcSTEP3: %d %lf %lf %lf %lf %lf %lf\n",i+1,FN->CPx[i].comp[0],FN->CPx[i].comp[1],FN->CPx[i].comp[2],FN->CPv[i].comp[0],FN->CPv[i].comp[1],FN->CPv[i].comp[2]);
    
  loop(i,nb) {
    loop(j,4) {
      int cpid = FN->Bz_list[i].CP[j]-1;
      loop(k,3) {
        FN->Bz_list[i].x[j*3+k] = FN->CPx[cpid].comp[k];
        FN->Bz_list[i].w[j*3+k] = FN->CPw[cpid].comp[k];
	  }
    }
    FN->Bz_list[i].update_bezier(1);
  }
  
  loop(i,np) loop(j,3) {
    gsl_matrix_set(FN->Psys, i,j, FN->CPx[i].comp[j]);
    gsl_matrix_set(FN->Wsys, i,j, FN->CPw[i].comp[j]);
  }
 
  //gsl_matrix_print(FN->Psys,"Psys:");
  //gsl_matrix_print(FN->Wsys,"Wsys:"); 
  //FN->assemble_matrices(); 

  //FN->compute_dH(0);
  
  //apply monotonicity constraints on bezier CPs
  //if(0) {
    //FN->compute_dH(0);
    //loop(i,nb) {
      //int cp[4] = { FN->Bz_list[i].CP[0]-1, FN->Bz_list[i].CP[1]-1, FN->Bz_list[i].CP[2]-1, FN->Bz_list[i].CP[3]-1 };
      ////printf("\n %d %d %d %d",cp[0],cp[1],cp[2],cp[3]);
      
      //double r01[3],r23[3],r12[3],r03[3],c01=0,c23=0,c12=0,l03=0;
    
      //loop(i,3) {
        //r01[i] = FN->CPx[cp[1]].comp[i] - FN->CPx[cp[0]].comp[i]; 
        //r23[i] = FN->CPx[cp[3]].comp[i] - FN->CPx[cp[2]].comp[i]; 
        //r12[i] = FN->CPx[cp[2]].comp[i] - FN->CPx[cp[1]].comp[i]; 
        //r03[i] = FN->CPx[cp[3]].comp[i] - FN->CPx[cp[0]].comp[i]; 
        
        //c01 += r01[i]*r03[i]; c23 += r23[i]*r03[i]; c12 += r12[i]*r03[i];
        //l03 += r03[i]*r03[i];
      //}
    
      //if(c01<0) { //CP1 has drifted to the "left" of CP0
        
        //printf("CP1 has drifted to the 'left' of CP0, adjusting .\n");
        //double vr=0;
        //loop(i,3) vr += FN->CPv[cp[1]].comp[i]*r03[i];
        
        //double fac = 2*fabs(c01)/l03, fac2;
        //if(vr<0) fac2 = 2*fabs(vr)/l03;
        
        //loop(i,3) {
          //FN->CPx[cp[1]].comp[i] += fac*r03[i];
          
          //if(vr<0) FN->CPv[cp[1]].comp[i] += fac2*r03[i];
        //}
      //}
      
      //if(c23<0) { //CP2 has drifted to the "right" of CP3
        
        //printf("CP2 has drifted to the 'right' of CP3, adjusting ...\n");
        //double vr=0;
        //loop(i,3) vr += FN->CPv[cp[2]].comp[i]*r03[i];
        
        //double fac = 2*fabs(c23)/l03, fac2;
        //if(vr>0) fac2 = 2*fabs(vr)/l03;
        
        //loop(i,3) {
          //FN->CPx[cp[2]].comp[i] -= fac*r03[i];
          
          //if(vr>0) FN->CPv[cp[2]].comp[i] -= fac2*r03[i];
        //}
      //}
      
      //if(c12<0) { //CP1 has drifted to the "right" of CP2
        
        //printf("CP1 has drifted to the 'right' of CP2, adjusting ...\n");
        //double vr1=0,vr2;
        //loop(i,3) {
          //vr1 += FN->CPv[cp[1]].comp[i]*r03[i];
          //vr2 += FN->CPv[cp[2]].comp[i]*r03[i];
        //}
        
        //double fac = fabs(c12)/l03, fac1,fac2;
        //if(vr1>0) fac1 = 2*fabs(vr1)/l03;
        //if(vr2<0) fac2 = 2*fabs(vr2)/l03;
        
        //loop(i,3) {
          //FN->CPx[cp[2]].comp[i] += fac*r03[i];
          //FN->CPx[cp[1]].comp[i] -= fac*r03[i];
          
          //if(vr1>0) FN->CPv[cp[1]].comp[i] -= fac1*r03[i];
          //if(vr2<0) FN->CPv[cp[2]].comp[i] += fac2*r03[i];
          
        //}
      //}
    //}
    
    //loop(i,np) loop(j,3) {
      //gsl_matrix_set(FN->Vsys, i,j, FN->CPv[i].comp[j]);
    //}
    
    //gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, FN->GM, FN->Vsys, 0.0, FN->Wsys);
  //}
  
  FN->compute_dH(1);
}

void apply_constraints(double q[], double p[], void *fn_data)
{
  fibnetwork *FN = (fibnetwork *) fn_data;

  int np = FN->np;	int nb = FN->nb;  int nf = FN->nf;
  
  loop(i,np) loop(j,3) if(!isnan(q[i*3+j])) FN->CPx[i].comp[j] = q[i*3+j];
  //loop(i,np) loop(j,3) if(!isnan(p[i*3+j])) FN->CPw[i].comp[j] = p[i*3+j]; //define momenta from momenta vector, define velocity later
  
  FN->compute_dH(0);
  
  //apply fixed constraints
  for(std::map<int,double>::iterator it = FN->fixed_dof.begin(); it != FN->fixed_dof.end(); it++) {
    int index = it->first-1;
    int i0 = index/3;
    int j0 = index%3;
    FN->CPx[i0].comp[j0] = it->second;
    FN->CPv[i0].comp[j0] = 0.0;
    //printf("FC APPLY: %d %d %d %lf %lf\n",i0,j0,index,it->second,FN->CPx[i0].comp[j0]);
  }
  
  //apply C1/C2
  if(1)
    loop(i,FN->C2_cps.size()) {
      int a = FN->C2_cps[i].comp[0];
      int b = FN->C2_cps[i].comp[1];
      int c = FN->C2_cps[i].comp[2];
      int d = FN->C2_cps[i].comp[3];
      int e = FN->C2_cps[i].comp[4];
      
      Vector3 db(FN->CPx[d-1].comp[0]-FN->CPx[b-1].comp[0], FN->CPx[d-1].comp[1]-FN->CPx[b-1].comp[1] , FN->CPx[d-1].comp[2]-FN->CPx[b-1].comp[2] );
      Vector3 eb(FN->CPx[e-1].comp[0]-FN->CPx[b-1].comp[0], FN->CPx[e-1].comp[1]-FN->CPx[b-1].comp[1] , FN->CPx[e-1].comp[2]-FN->CPx[b-1].comp[2] );
      Vector3 ad(FN->CPx[a-1].comp[0]-FN->CPx[d-1].comp[0], FN->CPx[a-1].comp[1]-FN->CPx[d-1].comp[1] , FN->CPx[a-1].comp[2]-FN->CPx[d-1].comp[2] );
      
      double m = sqrt(cross_mag(db,eb)/cross_mag(db,ad));
      double x = 1.0/(1.0+m);
      
      //EDIT to avoid singularities from control point overlap
      if(x<0.05) x = 0.05; if(x>0.95) x = 0.95;
      
      //x = 0.5;
      
      //printf("%d %d %d %d %d %d %lf %lf\n",i,a,b,c,d,e,m,x);
      loop(j,3) FN->CPx[c-1].comp[j] = x*FN->CPx[d-1].comp[j] + (1.0-x)*FN->CPx[b-1].comp[j];
    }
  
  //loop(i,np) printf("funcSTEP3: %d %lf %lf %lf %lf %lf %lf\n",i+1,FN->CPx[i].comp[0],FN->CPx[i].comp[1],FN->CPx[i].comp[2],FN->CPv[i].comp[0],FN->CPv[i].comp[1],FN->CPv[i].comp[2]);
  //printf("P Before \n"); loop(i,np) printf("%lf %lf %lf\n",p[i*3],p[i*3+1],p[i*3+2]);
  
  if(1) {
    
    //FN->compute_dH(0);
    
    loop(i,nb) {
      int cp[4] = { FN->Bz_list[i].CP[0]-1, FN->Bz_list[i].CP[1]-1, FN->Bz_list[i].CP[2]-1, FN->Bz_list[i].CP[3]-1 };
      //printf("\n %d %d %d %d",cp[0],cp[1],cp[2],cp[3]);
    
      Vector3 c03(FN->CPx[cp[3]].comp[0] - FN->CPx[cp[0]].comp[0],
                  FN->CPx[cp[3]].comp[1] - FN->CPx[cp[0]].comp[1],
                  FN->CPx[cp[3]].comp[2] - FN->CPx[cp[0]].comp[2]);
      
      Vector3 c01(FN->CPx[cp[1]].comp[0] - FN->CPx[cp[0]].comp[0],
                  FN->CPx[cp[1]].comp[1] - FN->CPx[cp[0]].comp[1],
                  FN->CPx[cp[1]].comp[2] - FN->CPx[cp[0]].comp[2]);
      
      Vector3 c23(FN->CPx[cp[3]].comp[0] - FN->CPx[cp[2]].comp[0],
                  FN->CPx[cp[3]].comp[1] - FN->CPx[cp[2]].comp[1],
                  FN->CPx[cp[3]].comp[2] - FN->CPx[cp[2]].comp[2]);
      
      double c03_1 = cross_mag(c01,c03)/(vec_mag(c01)*vec_mag(c03));
      double c03_2 = cross_mag(c23,c03)/(vec_mag(c23)*vec_mag(c03));
      
      printf("nbE %16.6f \t Ltest %16.6f %16.6f\n",0.0,c03_1,c03_2);
      
      //double r01[3],r23[3],r12[3],r03[3],c01=0,c23=0,c12=0,l03=0;
        
      //if(c03_1<BEZ_EPS) 
      if(1) {
        
        double r03[3],r13[3],c13=0,l03=0,vr=0;
        
        loop(j,3) {
          r13[j] = FN->CPx[cp[3]].comp[j] - FN->CPx[cp[1]].comp[j]; 
          r03[j] = FN->CPx[cp[3]].comp[j] - FN->CPx[cp[0]].comp[j]; 
          
          vr += FN->CPv[cp[1]].comp[j]*r03[j];
          c13 += r13[j]*r03[j]; l03 += r03[j]*r03[j];
        }
        
        if(c13<0) {//CP1 has drifted to the right of CP3, bringing back
          
          double fac = 2*fabs(c13)/l03, fac2 = 2*fabs(vr)/l03;
          
          loop(j,3) {            
            //FN->CPx[cp[1]].comp[j] -= fac*r03[j];
            if(vr>0) FN->CPv[cp[1]].comp[j] -= fac2*r03[j];
            
            //FN->CPx[cp[1]].comp[j] = 0.33*FN->CPx[cp[3]].comp[j] + 0.67*FN->CPx[cp[0]].comp[j];
            //if(j==0) FN->CPv[cp[1]].comp[j] *= -1.0;
        
          }
        }
      }
      
      //if(c03_2<BEZ_EPS) 
      if(1) {
        
        double r03[3],r02[3],c02=0,l03=0,vr=0;
        
        loop(j,3) {
          r02[j] = FN->CPx[cp[2]].comp[j] - FN->CPx[cp[0]].comp[j]; 
          r03[j] = FN->CPx[cp[3]].comp[j] - FN->CPx[cp[0]].comp[j]; 
          
          vr += FN->CPv[cp[2]].comp[j]*r03[j];
          c02 += r02[j]*r03[j]; l03 += r03[j]*r03[j];
        }
        
        if(c02<0) {//CP2 has drifted to the left of CP0, bringing back
          
          double fac = 2*fabs(c02)/l03, fac2 = 2*fabs(vr)/l03;
          
          loop(j,3) {            
            //FN->CPx[cp[2]].comp[j] += fac*r03[j];
            if(vr<0) FN->CPv[cp[2]].comp[j] += fac2*r03[j];
            
            //FN->CPx[cp[1]].comp[j] = 0.33*FN->CPx[cp[3]].comp[j] + 0.67*FN->CPx[cp[0]].comp[j];
            //if(j==0) FN->CPv[cp[1]].comp[j] *= -1.0;
        
          }
        }
      }

      //if(c03_2<BEZ_EPS) loop(j,3) {
          ////FN->CPx[cp[2]].comp[j] = 0.33*FN->CPx[cp[0]].comp[j] + 0.67*FN->CPx[cp[3]].comp[j];
          //if(j==0) FN->CPv[cp[2]].comp[j] *= -1.0;
      //}
      
    }
    
  }
  
  if(1) {
    
    //printf("V After \n");
    //loop(i,np) printf("%lf %lf %lf\n",FN->CPv[i].comp[0],FN->CPv[i].comp[1],FN->CPv[i].comp[2]);
  
    loop(i,np) loop(j,3) {
      gsl_matrix_set(FN->Vsys, i,j, FN->CPv[i].comp[j]);
    }
    
    //gsl_matrix_print(FN->Wsys,"Wsys before:");
    gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, FN->GM, FN->Vsys, 0.0, FN->Wsys);
    
    //gsl_matrix_print(FN->GM,"ac GM:");
    //gsl_matrix_print(FN->GMinv,"ac GMinv:");
    
    //gsl_matrix_print(FN->Wsys,"Wsys after:");
    
    loop(i,np) loop(j,3) {
      FN->CPw[i].comp[j] = gsl_matrix_get(FN->Wsys, i,j);
    }
  
    loop(i,np) loop(j,3) if(!isnan( FN->CPw[i].comp[j] )) p[i*3+j] = FN->CPw[i].comp[j];
  }
  
  if(0) {
    
    FN->compute_dH(0);
    
    //printf("V Before \n");
    //loop(i,np) printf("%lf %lf %lf\n",FN->CPv[i].comp[0],FN->CPv[i].comp[1],FN->CPv[i].comp[2]);
    
    loop(i,nb) {
      int cp[4] = { FN->Bz_list[i].CP[0]-1, FN->Bz_list[i].CP[1]-1, FN->Bz_list[i].CP[2]-1, FN->Bz_list[i].CP[3]-1 };
      //printf("\n %d %d %d %d",cp[0],cp[1],cp[2],cp[3]);
      
      double r01[3],r23[3],r12[3],r03[3],c01=0,c23=0,c12=0,l03=0;
    
      loop(i,3) {
        r01[i] = FN->CPx[cp[1]].comp[i] - FN->CPx[cp[0]].comp[i]; 
        r23[i] = FN->CPx[cp[3]].comp[i] - FN->CPx[cp[2]].comp[i]; 
        r12[i] = FN->CPx[cp[2]].comp[i] - FN->CPx[cp[1]].comp[i]; 
        r03[i] = FN->CPx[cp[3]].comp[i] - FN->CPx[cp[0]].comp[i]; 
        
        c01 += r01[i]*r03[i]; c23 += r23[i]*r03[i]; c12 += r12[i]*r03[i];
        l03 += r03[i]*r03[i];
      }
    
      if(1) if(c01<0) { //CP1 has drifted to the "left" of CP0
        
        printf("CP1 has drifted to the 'left' of CP0, adjusting ... %lf\n", c01);
        double vr=0;
        loop(i,3) vr += FN->CPv[cp[1]].comp[i]*r03[i];
        
        double fac = 2*fabs(c01)/l03, fac2;
        if(vr<0) fac2 = 2*fabs(vr)/l03;
        
        loop(i,3) {
          FN->CPx[cp[1]].comp[i] += fac*r03[i];
          
          if(vr<0) FN->CPv[cp[1]].comp[i] += fac2*r03[i];
        }
      }
      
      if(1) if(c23<0) { //CP2 has drifted to the "right" of CP3
        
        printf("CP2 has drifted to the 'right' of CP3, adjusting ... %lf\n", c23);
        double vr=0;
        loop(i,3) vr += FN->CPv[cp[2]].comp[i]*r03[i];
        
        double fac = 2*fabs(c23)/l03, fac2;
        if(vr>0) fac2 = 2*fabs(vr)/l03;
        
        loop(i,3) {
          FN->CPx[cp[2]].comp[i] -= fac*r03[i];
          
          if(vr>0) FN->CPv[cp[2]].comp[i] -= fac2*r03[i];
        }
      }
      
      if(1) if(c12<0) { //CP1 has drifted to the "right" of CP2
        
        printf("CP1 has drifted to the 'right' of CP2, adjusting .. %lf\n", c12);
        double vr1=0,vr2=0;
        loop(i,3) {
          vr1 += FN->CPv[cp[1]].comp[i]*r03[i];
          vr2 += FN->CPv[cp[2]].comp[i]*r03[i];
        }
        
        double fac = fabs(c12)/l03, fac1,fac2;
        if(vr1>0) fac1 = 2*fabs(vr1)/l03;
        if(vr2<0) fac2 = 2*fabs(vr2)/l03;
        
        loop(i,3) {
          FN->CPx[cp[2]].comp[i] += fac*r03[i];
          FN->CPx[cp[1]].comp[i] -= fac*r03[i];
          
          //printf("v1 before: %lf %lf %lf\n",FN->CPv[cp[1]].comp[0],FN->CPv[cp[1]].comp[1],FN->CPv[cp[1]].comp[2]);
          //printf("v2 before: %lf %lf %lf\n",FN->CPv[cp[2]].comp[0],FN->CPv[cp[2]].comp[1],FN->CPv[cp[2]].comp[2]);
          
          //if(vr1>0) FN->CPv[cp[1]].comp[i] -= fac1*r03[i];
          //if(vr2<0) FN->CPv[cp[2]].comp[i] += fac2*r03[i];
          
          if(vr1>0) FN->CPv[cp[1]].comp[i] = 0.0;
          if(vr2<0) FN->CPv[cp[2]].comp[i] = 0.0;
          
          //printf("v1 after: %lf %lf %lf\n",FN->CPv[cp[1]].comp[0],FN->CPv[cp[1]].comp[1],FN->CPv[cp[1]].comp[2]);
          //printf("v2 after: %lf %lf %lf\n",FN->CPv[cp[2]].comp[0],FN->CPv[cp[2]].comp[1],FN->CPv[cp[2]].comp[2]);

          
        }
      }
    }
    
    if(1) {
      
      //printf("V After \n");
      //loop(i,np) printf("%lf %lf %lf\n",FN->CPv[i].comp[0],FN->CPv[i].comp[1],FN->CPv[i].comp[2]);
    
      loop(i,np) loop(j,3) {
        gsl_matrix_set(FN->Vsys, i,j, FN->CPv[i].comp[j]);
      }
      
      //gsl_matrix_print(FN->Wsys,"Wsys before:");
      gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, FN->GM, FN->Vsys, 0.0, FN->Wsys);
      
      //gsl_matrix_print(FN->GM,"ac GM:");
      //gsl_matrix_print(FN->GMinv,"ac GMinv:");
      
      
      
      //gsl_matrix_print(FN->Wsys,"Wsys after:");
      
      loop(i,np) loop(j,3) {
        FN->CPw[i].comp[j] = gsl_matrix_get(FN->Wsys, i,j);
      }
    
      loop(i,np) loop(j,3) if(!isnan( FN->CPw[i].comp[j] )) p[i*3+j] = FN->CPw[i].comp[j];
    }
  }
  
  //printf("P After \n"); loop(i,np) printf("%lf %lf %lf\n",p[i*3],p[i*3+1],p[i*3+2]);
  
  loop(i,np) loop(j,3) if(!isnan( FN->CPx[i].comp[j] )) q[i*3+j] = FN->CPx[i].comp[j];
  
  
  //loop(i,np) loop(j,3) if(!isnan( FN->CPw[i].comp[j] )) p[i*3+j] = FN->CPw[i].comp[j]; //define momenta from momenta vector, define velocity later
 
  
  //loop(i,nb) {
    //loop(j,4) {
      //int cpid = FN->Bz_list[i].CP[j]-1;
      //loop(k,3) {
        //FN->Bz_list[i].x[j*3+k] = FN->CPx[cpid].comp[k];
        //FN->Bz_list[i].w[j*3+k] = FN->CPw[cpid].comp[k];
	  //}
    //}
    //FN->Bz_list[i].update_bezier(1);
  //}
  
  //loop(i,np) loop(j,3) {
    //gsl_matrix_set(FN->Psys, i,j, FN->CPx[i].comp[j]);
    //gsl_matrix_set(FN->Wsys, i,j, FN->CPw[i].comp[j]);
  //}
  
  //FN->assemble_matrices(); 

  //FN->compute_dH();
}

void tao_integrate(fibnetwork &fn, int nsteps, double tstep)
{
  size_t np = fn.np;
	
	fn.Psys = gsl_matrix_alloc (np, 3);
	fn.Qsys = gsl_matrix_alloc (np, 3);
	fn.Rsys = gsl_matrix_alloc (np, 3);
	fn.Fsys = gsl_matrix_alloc (np, 3);
	fn.Wsys = gsl_matrix_alloc (np, 3);
  fn.Vsys = gsl_matrix_alloc (np, 3);
  fn.qH = gsl_matrix_alloc (np, 3);
  fn.pH = gsl_matrix_alloc (np, 3);
  
  
	fn.GM = gsl_matrix_alloc (np, np);
	fn.GMinv = gsl_matrix_alloc (np, np);
	fn.GMdot = gsl_matrix_alloc (np, np);
	fn.GN = gsl_matrix_alloc (np, np);
  
  //////
  //fn2.Psys = gsl_matrix_alloc (np, 3);
	//fn2.Qsys = gsl_matrix_alloc (np, 3);
	//fn2.Rsys = gsl_matrix_alloc (np, 3);
	//fn2.Fsys = gsl_matrix_alloc (np, 3);
	//fn2.Wsys = gsl_matrix_alloc (np, 3);
  //fn2.qH = gsl_matrix_alloc (np, 3);
  //fn2.pH = gsl_matrix_alloc (np, 3);
  
	//fn2.GM = gsl_matrix_alloc (np, np);
  //fn2.GMinv = gsl_matrix_alloc (np, np);
	//fn2.GMdot = gsl_matrix_alloc (np, np);
	//fn2.GN = gsl_matrix_alloc (np, np);
  //////
  
  double *X = new double[np*3];
  double *Y = new double[np*3];
  double *Q = new double[np*3];
  double *P = new double[np*3];
  
	loop(i,np) loop(j,3) {
    
    Q[i*3+j] = fn.CPx[i].comp[j];
    Y[i*3+j] = fn.CPw[i].comp[j];
    
    X[i*3+j] = fn.CPx[i].comp[j];
    P[i*3+j] = fn.CPw[i].comp[j];

  }

  int ac = 1;
  //print_QPXY(Q,P,X,Y,np);

  for (int s = 1; s <= nsteps; s++) {
    
    //printf("\nSTART TAO STEP:\n");
    //dH(Q,P,&fn,0);
    //fn.totE = fn.axialE + fn.bendE + fn.cohE + fn.kinE;
    //printf("\nPRE: %d %lf %lf %lf %lf %lf\n", s, fn.axialE, fn.bendE, fn.cohE, fn.kinE, fn.totE);
    
    //if(0) { //do not apply constraint here as matrices are not defined at this stage
      //apply_constraints(Q,P,&fn);
    	//apply_constraints(X,Y,&fn);
    //}
    //STEP 1
    //printf("\nSTEP 1:\n");
    dH(Q,Y,&fn,0);
    loop(i,np) loop(j,3) {
      P[i*3+j] -= 0.5*tstep*gsl_matrix_get(fn.qH,i,j);
      X[i*3+j] += 0.5*tstep*gsl_matrix_get(fn.pH,i,j);
    }
    //apply_constraints(X,Y,&fn);
    //print_QPXY(Q,P,X,Y,np);
 
    //STEP 2
    //printf("\nSTEP 2:\n");
    dH(X,P,&fn,0);
    loop(i,np) loop(j,3) {
      Q[i*3+j] += 0.5*tstep*gsl_matrix_get(fn.pH,i,j);
      Y[i*3+j] -= 0.5*tstep*gsl_matrix_get(fn.qH,i,j);
    }
    //apply_constraints(Q,P,&fn);
    //print_QPXY(Q,P,X,Y,np);

    if(ac) {
        apply_constraints(Q,P,&fn);
        apply_constraints(X,Y,&fn);
    }
    
    //STEP 3
    //printf("\nSTEP 3:\n");
    loop(i,np) loop(j,3) {
      
      double a = Q[i*3+j]+X[i*3+j];
      double b = Q[i*3+j]-X[i*3+j];
      
      double c = P[i*3+j]+Y[i*3+j];
      double d = P[i*3+j]-Y[i*3+j];
      
      Q[i*3+j] = 0.5*( a + b*cos(2*tao_omega*tstep) + d*sin(2*tao_omega*tstep) );
      P[i*3+j] = 0.5*( c + d*cos(2*tao_omega*tstep) - b*sin(2*tao_omega*tstep) );
      X[i*3+j] = 0.5*( a - b*cos(2*tao_omega*tstep) - d*sin(2*tao_omega*tstep) );
      Y[i*3+j] = 0.5*( c - d*cos(2*tao_omega*tstep) + b*sin(2*tao_omega*tstep) );
    }
    //print_QPXY(Q,P,X,Y,np);
    //apply_constraints(Q,P,&fn);
    //apply_constraints(X,Y,&fn);
    
    if(ac) {
        apply_constraints(Q,P,&fn);
        apply_constraints(X,Y,&fn);
    } 
    //STEP 4
    //printf("\nSTEP 4:\n");
    dH(X,P,&fn,0);
    loop(i,np) loop(j,3) {
      Q[i*3+j] += 0.5*tstep*gsl_matrix_get(fn.pH,i,j);
      Y[i*3+j] -= 0.5*tstep*gsl_matrix_get(fn.qH,i,j);
    }
    //apply_constraints(Q,P,&fn);
    //print_QPXY(Q,P,X,Y,np);
    if(ac) {
        apply_constraints(Q,P,&fn);
        apply_constraints(X,Y,&fn);
    }


    //STEP 5
    //printf("\nSTEP 5:\n");
    dH(Q,Y,&fn,0);
    loop(i,np) loop(j,3) {
      P[i*3+j] -= 0.5*tstep*gsl_matrix_get(fn.qH,i,j);
      X[i*3+j] += 0.5*tstep*gsl_matrix_get(fn.pH,i,j);
    }
    //apply_constraints(X,Y,&fn);
    //print_QPXY(Q,P,X,Y,np);

    //print stats
    //printf("\nSTEP STATS:\n");
    dH(Q,P,&fn,0);
    dH(X,Y,&fn,0);
    
    if(ac) {
        apply_constraints(Q,P,&fn);
        apply_constraints(X,Y,&fn);
    }
    
    //printf("\nSTEP 5.1:\n");
    dH(Q,P,&fn,0);
   
    fn.tstep = s;	 
    fn.printlammps_cps((char*) "Output_cps.lammpstrj",(char*) "a");
    fn.printlammps((char*) "Output.lammpstrj",(char*) "a",NBEZ);
    
    fn.totE = fn.axialE + fn.bendE + fn.cohE + fn.kinE;
    printf("tao_int %d %lf %lf %lf %lf %lf : %.12f %.12f %.12f\n", s, fn.axialE, fn.bendE, fn.cohE, fn.kinE, fn.totE, fn.momT[0], fn.momT[1], fn.momT[2]);
  }
  
  delete[] X;
  delete[] Y;
  delete[] Q;
  delete[] P;
  
  gsl_matrix_free(fn.Psys);
  gsl_matrix_free(fn.Qsys);
  gsl_matrix_free(fn.Rsys);
  gsl_matrix_free(fn.Fsys);
  gsl_matrix_free(fn.Wsys);
  gsl_matrix_free(fn.Vsys);
  gsl_matrix_free(fn.qH);
  gsl_matrix_free(fn.pH);
  
  gsl_matrix_free(fn.GM);
  gsl_matrix_free(fn.GMinv);
  gsl_matrix_free(fn.GMdot);
  gsl_matrix_free(fn.GN);
  
  //gsl_matrix_free(fn2.Psys);
  //gsl_matrix_free(fn2.Qsys);
  //gsl_matrix_free(fn2.Rsys);
  //gsl_matrix_free(fn2.Fsys);
  //gsl_matrix_free(fn2.Wsys);
  //gsl_matrix_free(fn2.qH);
  //gsl_matrix_free(fn2.pH);
  
  //gsl_matrix_free(fn2.GM);
  //gsl_matrix_free(fn2.GMinv);
  //gsl_matrix_free(fn2.GMdot);
  //gsl_matrix_free(fn2.GN);
  
}
  
