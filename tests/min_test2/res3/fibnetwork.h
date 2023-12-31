class fibnetwork {
  
  public:
    int np,nb,nf,nt,ncc,tstep;
    double bounds[3][2];
    double axialE, bendE, cohE, kinE, totE, time, momT[3];
    
    std::map<int,double> fixed_dof;
    std::vector<iVector3> C1_cps;
    std::vector<iVector5> C2_cps;
    std::vector<Vector3> CPx,CPv,CPf,CPa,CPw,CPdof;
    std::vector<cBezier> Bz_list;
    std::vector<std::vector<int>> fibers, CP_bid; 
    std::vector<int> CP_fid; 
    
    gsl_matrix *GM, *GMdot, *GN, *GMinv;
    gsl_matrix *Psys, *Qsys, *Rsys, *Fsys, *Wsys, *qH, *pH, *Vsys;
    
    void readfile(char *);
    void printdata();
    void printlammps(char *,char *, int);
    void printlammps_cps(char *,char *);
    void compute_ef();
    void assemble_matrices();
    void constraint_energy();
    void minimize();
    void computeK(gsl_matrix *);
    
    double dMdQ(int, int, int, int);
    void compute_dH(int);
    //void update_sys(double,gsl_matrix *,gsl_matrix *, int);
    //void integrate_runge_kutta_4(int, double);  
};

void fibnetwork::readfile(char *filename)
{
  std::ifstream infile(filename);
  
  std::string word,line;
  loop(i,2) std::getline(infile, line);
  
  sss ss0;
  std::getline(infile, line); ss0 << line; ss0 >> np; ss0.str("");
  std::getline(infile, line); ss0 << line; ss0 >> nb; ss0.str("");
  std::getline(infile, line); ss0 << line; ss0 >> nf; ss0.str("");
  //std::getline(infile, line); ss0 << line; ss0 >> ncc; ss0.str("");
  std::getline(infile, line);
  std::getline(infile, line); ss0 << line; ss0 >> nt; ss0.str("");
  std::getline(infile, line);
  
  loop(i,3) {
    sss ss;
    std::getline(infile, line); ss << line; 
    ss >> bounds[i][0] >> bounds[i][1];    
  }
  loop(i,3) std::getline(infile, line);
  
  loop(i,BTYPES) {
    sss ss;
    std::getline(infile, line); ss << line;
    int ctmp;  ss >> ctmp >> EA[ctmp-1] >> EI[ctmp-1] >> rho0[ctmp-1]; ss.str("");
    //printf("%s %d %lf %lf\n",line.c_str(),ctmp,EA[ctmp-1],EI[ctmp-1]);
  }
  
  loop(i,3) std::getline(infile, line);

  std::srand(std::time(0)); 
  loop(i,np) {
    sss ss;
    double x[3],v[3],tmp0[3]; int tmp,fdof[3];
    std::getline(infile, line); ss << line; 
    ss >> tmp >> tmp0[0] >> tmp0[1] >> tmp0[2] >> x[0] >> x[1] >> x[2] >> v[0] >> v[1] >> v[2] >> fdof[0] >> fdof[1] >> fdof[2];
    
    loop(j,3) if(fdof[j]==0) fixed_dof.insert(std::make_pair((tmp-1)*3+j+1, x[j] ));
    
    double rn = RTOL*((double) std::rand() / (RAND_MAX));
    
    CPx.push_back(Vector3(x[0]+rn,x[1]+rn,x[2]+rn));
    CPw.push_back(Vector3(v[0],v[1],v[2]));
    CPdof.push_back(Vector3(fdof[0],fdof[1],fdof[2]));
    CPv.push_back(Vector3(0,0,0));
    CPf.push_back(Vector3(0,0,0));
    CPa.push_back(Vector3(0,0,0));
    CP_fid.push_back(0);
    
    std::vector<int> bid_rank;
    CP_bid.push_back(bid_rank);
  }
  loop(i,3) std::getline(infile, line);
  
  loop(i,nb) {
    sss ss;
    cBezier bzi; int tmp[6]; double length0; 
    std::getline(infile, line); ss << line; 
    ss >> tmp[0] >> tmp[1] >> tmp[2] >> tmp[3] >> tmp[4] >> tmp[5] >> length0;
    
    bzi.type = tmp[1]; bzi.id = tmp[0];
    loop(j,4) {
      bzi.CP[j] = tmp[j+2]; 
      int cpid = bzi.CP[j]-1;
      loop(k,3) bzi.x[k+j*3] = CPx[cpid].comp[k];
    }
    bzi.length0 = length0;
    //printf("Initialize bezier: %d ",i+1);
    bzi.initialize_bezier();
    
    Bz_list.push_back(bzi);
    
    loop(j,4) {
      CP_bid[(tmp[j+2]-1)].push_back(tmp[0]); //CP_bid[(tmp[j+2]-1)].push_back(j+1);
    }
  }
  loop(i,3) std::getline(infile, line);
  
  loop(i,nf) {
    sss ss; int tmp,fid;
    std::vector<int> fbz;
    std::getline(infile, line); ss << line; ss >> fid;
    while (ss >> tmp) fbz.push_back(tmp);
    
    loop(j,fbz.size()) {
      Bz_list[fbz[j]-1].fiber = fid;
      
      //printf("Fiber:%d CPs:%d %d %d %d\n",fid,Bz_list[fbz[j]-1].CP[0],Bz_list[fbz[j]-1].CP[1],Bz_list[fbz[j]-1].CP[2],Bz_list[fbz[j]-1].CP[3]);
      
      loop(k,4) {
        int tmpk = Bz_list[fbz[j]-1].CP[k];
        CP_fid[tmpk-1] = fid;
      }
      
      if(j==0) {Bz_list[fbz[j]-1].leftB = -1; Bz_list[fbz[j]-1].rightB = fbz[j+1];}
      else if(j== fbz.size()-1) {Bz_list[fbz[j]-1].leftB = fbz[j-1]; Bz_list[fbz[j]-1].rightB = -1;}
      else {Bz_list[fbz[j]-1].leftB = fbz[j-1]; Bz_list[fbz[j]-1].rightB = fbz[j+1];}
    } 
    
    ncc = 0;
    loop(j,fbz.size()) {
      if(j==0) continue;
      int lb = Bz_list[fbz[j]-1].leftB;
      int curr = fbz[j];
            
      C1_cps.push_back(iVector3(Bz_list[lb-1].CP[2], Bz_list[lb-1].CP[3], Bz_list[curr-1].CP[1]));
      C2_cps.push_back(iVector5(Bz_list[lb-1].CP[1], Bz_list[lb-1].CP[2], Bz_list[lb-1].CP[3], Bz_list[curr-1].CP[1], Bz_list[curr-1].CP[2]));
      ncc++;
   }
    
    fibers.push_back(fbz);
  }
  
  CPv.resize(CPx.size()); CPf.resize(CPx.size());
  CPw.resize(CPx.size());
  
  
  infile.close();
}

void fibnetwork::printdata()
{
  printf("Control points \n");
  loop(i,np) {
    printf("%d %lf %lf %lf \n",i+1, CPx[i].comp[0], CPx[i].comp[1], CPx[i].comp[2]);
  }
  
  printf("\nCubic Beziers \n");
  loop(i,nb) {
    printf("%d %d %d %d %d %d (%d,%d)\n",Bz_list[i].id, Bz_list[i].type, Bz_list[i].CP[0], Bz_list[i].CP[1], Bz_list[i].CP[2], Bz_list[i].CP[3], Bz_list[i].leftB, Bz_list[i].rightB);
  }
  
  printf("\nFibers \n");
  loop(i,nf) {
    std::cout << i+1 << ' ';
    for (auto k: fibers[i]) std::cout << k << ' ';
    std::cout << std::endl ;
  }
}

void fibnetwork::printlammps_cps(char *filename, char *mode)
{
  FILE *fp = fopen(filename, mode);
  
  fprintf(fp, "ITEM: TIMESTEP\n%d\nITEM: NUMBER OF ATOMS\n%d\n",tstep,np);
  fprintf(fp, "ITEM: BOX BOUNDS pp pp pp\n%2.6f %2.6f \n%2.6f %2.6f \n%2.6f %2.6f \n", 
      bounds[0][0], bounds[0][1], bounds[1][0], bounds[1][1], bounds[2][0], bounds[2][1]);
  fprintf(fp, "ITEM: ATOMS id mol type q xu yu zu px py pz\n");
  
  loop(i,np) {
    int bezid = CP_bid[i][0];
    fprintf(fp, "%d %d 1 %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.0f %.0f %.0f\n", i+1, CP_fid[i], Bz_list[bezid-1].length0, CPx[i].comp[0], CPx[i].comp[1], CPx[i].comp[2], CPw[i].comp[0], CPw[i].comp[1], CPw[i].comp[2], CPdof[i].comp[0], CPdof[i].comp[1], CPdof[i].comp[2] );
  }
  
  fclose(fp);
}

void fibnetwork::printlammps(char *filename, char *mode, int n)
{
  //int num_coord = 3*(n*nb+1);
  
  std::vector<Vector4> points;
  //points.push_back(Vector3(Bz_list[fibers[0][0]-1].x[0],Bz_list[fibers[0][0]-1].x[1],Bz_list[fibers[0][0]-1].x[2]));
  
  loop(i,nf) {
    points.push_back(Vector4(Bz_list[fibers[i][0]-1].fiber,Bz_list[fibers[i][0]-1].x[0],Bz_list[fibers[i][0]-1].x[1],Bz_list[fibers[i][0]-1].x[2]));
  
    for (auto k: fibers[i])
    {
      Bz_list[k-1].interpolateL(points,n);
    }
  }  
  
  FILE *fp = fopen(filename, mode);
  
  fprintf(fp, "ITEM: TIMESTEP\n%d\nITEM: NUMBER OF ATOMS\n%ld\n",tstep,points.size());
  fprintf(fp, "ITEM: BOX BOUNDS pp pp pp\n%2.6f %2.6f \n%2.6f %2.6f \n%2.6f %2.6f \n", 
      bounds[0][0], bounds[0][1], bounds[1][0], bounds[1][1], bounds[2][0], bounds[2][1]);
  fprintf(fp, "ITEM: ATOMS id mol type q xu yu zu\n");
  
  loop(i,points.size()) {
    fprintf(fp, "%d %.0f 1 0.00 %.6f %.6f %.6f\n", i+1, points[i].comp[0], points[i].comp[1], points[i].comp[2], points[i].comp[3]);
  }
  
  fclose(fp);
}

void fibnetwork::compute_ef()
{
  axialE = 0.0; bendE = 0.0; cohE = 0.0; kinE = 0.0;
  loop(i,np) loop(j,3) CPf[i].comp[j] = 0.0;
  
  loop(j,3) momT[j] = 0;
  
  loop(i,nf)
    for (auto k: fibers[i])
      loop(b,4) loop(j,3) Bz_list[k-1].f[b*3+j] = 0.0;

  //loop(i,nb) loop2(j,i,nb)
      //if (Bz_list[i].fiber != Bz_list[j].fiber)
        //cohE += inter_energy(&Bz_list[i],&Bz_list[j]);
  
  loop(i,nf) {
    for (auto k: fibers[i])
    {
      //loop(b,4) loop(j,3) Bz_list[k-1].f[b*3+j] = 0.0;
      
      Bz_list[k-1].axial_ef();    
      //printf("Forces: \n"); loop(i,4) printf("%.4e %.4e %.4e \n", Bz_list[k-1].f[i*3],Bz_list[k-1].f[i*3+1],Bz_list[k-1].f[i*3+2]);
      Bz_list[k-1].bending_ef();
      
      Bz_list[k-1].penalty_ef();
      
      
      Bz_list[k-1].momentum();
      Bz_list[k-1].kin_energy();
      //printf("Forces: \n"); loop(i,4) printf("%.4e %.4e %.4e \n", Bz_list[k-1].f[i*3],Bz_list[k-1].f[i*3+1],Bz_list[k-1].f[i*3+2]);
      //printf("%d %d %lf %lf\n",i+1,k,Bz_list[k-1].axialE, Bz_list[k-1].bendE);
      
      //loop(b,4) {
        //int id = Bz_list[k-1].CP[b];
        //loop(j,3) CPf[id-1].comp[j] += Bz_list[k-1].f[b*3+j];
      //}
      
      //printf("Bezier %d Axial Energy: %lf, %lf %lf\n",k,Bz_list[k-1].axialE,Bz_list[k-1].length,Bz_list[k-1].length0);
      
      axialE += Bz_list[k-1].axialE;
      bendE += Bz_list[k-1].bendE;
      kinE += Bz_list[k-1].kinE;
      //cohE += Bz_list[k-1].cohE;
      
      loop(j,3) momT[j] += Bz_list[k-1].mom[j];
      
      //loop(b,4) 
        //printf("Bezier: %d %lf %lf %lf\n",k,Bz_list[k-1].f[b*3],Bz_list[k-1].f[b*3+1],Bz_list[k-1].f[b*3+2]);
    }
  }  

  loop(i,nb) loop2(j,i,nb)
      if (Bz_list[i].fiber != Bz_list[j].fiber)
        cohE += inter_energy(&Bz_list[i],&Bz_list[j]); 
  
  //loop(i,np) printf("compute_ef,CPf0: %d %lf %lf %lf\n",i+1,CPf[i].comp[0],CPf[i].comp[1],CPf[i].comp[2]);


  loop(i,nf)
    for (auto k: fibers[i])
      loop(b,4) {
        int id = Bz_list[k-1].CP[b];
          loop(j,3) if(!isnan(Bz_list[k-1].f[b*3+j])) CPf[id-1].comp[j] += Bz_list[k-1].f[b*3+j];
      }
      
  //loop(i,np) printf("compute_ef,CPf3: %d %lf %lf %lf\n",i+1,CPf[i].comp[0],CPf[i].comp[1],CPf[i].comp[2]);
  
  
  //totE = axialE + bendE + cohE + kinE;
  //constraint_energy();  
  //printf("Energies: %lf %lf %lf %lf %lf\n", axialE, bendE, cohE, kinE, totE);
  
}

/////////////////MINIMIZE//////////////////////////

double VdV(unsigned n, const double *input, double *grad, void *fn_data)
{
  fibnetwork *FN = (fibnetwork *) fn_data;
    
  loop(i,FN->np) loop(j,3) FN->CPx[i].comp[j] = input[i*3+j];
    
  loop(i,FN->nb) {
    
    loop(j,4) {
      
      int cpid = FN->Bz_list[i].CP[j]-1;
      loop(k,3) 
        FN->Bz_list[i].x[j*3+k] = FN->CPx[cpid].comp[k];
    }
    FN->Bz_list[i].update_bezier(0);
  }
  
  FN->compute_ef();
      
  if (grad) {
    loop(i,FN->np) loop(j,3) grad[i*3+j] = -FN->CPf[i].comp[j] ;
  }
  
  //FN->printlammps_cps((char*) "Output_cps.lammpstrj",(char*) "a");
  //FN->printlammps((char*) "Output.lammpstrj",(char*) "a",20);
  printf("Min Energies: %lf %lf %lf %lf %lf\n", FN->axialE, FN->bendE, FN->cohE, FN->kinE, FN->totE); 
  return (FN->axialE + FN->bendE + FN->cohE);
}

void fixedconstraints(unsigned m, double *result, unsigned n, const double* input, double* grad, void* fn_data)
{
  fibnetwork *FN = (fibnetwork *) fn_data;
  
  if (grad) loop(i,m*n) grad[i] = 0.0;
  
  int ci = 0;
  for(std::map<int,double>::iterator it = FN->fixed_dof.begin(); it != FN->fixed_dof.end(); it++) {
    result[ci] = input[it->first-1] - it->second;
    if (grad) grad[ci*n+(it->first-1)] = 1.0;
    ci++;
  }
}

void C1_continuity(unsigned m, double *result, unsigned n, const double* input, double* grad, void* fn_data)
{
  fibnetwork *FN = (fibnetwork *) fn_data;
  
  if (grad) loop(i,m*n) grad[i] = 0.0;
  
  loop(i,FN->C1_cps.size()) {
    int a = FN->C1_cps[i].comp[0];
    int b = FN->C1_cps[i].comp[1];
    int c = FN->C1_cps[i].comp[2];
    
    Vector3 p,q;
    
    loop(j,3) {
      p.comp[j] = input[(b-1)*3+j] - input[(a-1)*3+j];
      q.comp[j] = input[(c-1)*3+j] - input[(b-1)*3+j];
    }

    double pp=0,qq=0,pq=0;
    loop(j,3) {
      pp += p.comp[j]*p.comp[j];
      qq += q.comp[j]*q.comp[j];
      pq += p.comp[j]*q.comp[j];
    }
    
    result[i] = pq-sqrt(pp*qq); 
    
    if(grad) {
      loop(j,3) {
        grad[i*n + (a-1)*3+j] =  -q.comp[j]  + sqrt(qq/pp)*p.comp[j];
        grad[i*n + (c-1)*3+j] =   p.comp[j]  - sqrt(pp/qq)*q.comp[j];
        grad[i*n + (b-1)*3+j] =   (q.comp[j] - p.comp[j]) + sqrt(pp/qq)*q.comp[j] - sqrt(qq/pp)*p.comp[j];
      }
    }
  }
}

void C2_continuity(unsigned m, double *result, unsigned n, const double* input, double* grad, void* fn_data)
{
  fibnetwork *FN = (fibnetwork *) fn_data;
  
  if (grad) loop(i,m*n) grad[i] = 0.0;
  
  loop(i,FN->C2_cps.size()) {
    int a = FN->C2_cps[i].comp[0];
    int b = FN->C2_cps[i].comp[1];
    int c = FN->C2_cps[i].comp[2];
    int d = FN->C2_cps[i].comp[3];
    int e = FN->C2_cps[i].comp[4];
    
    Vector3 p,q,A,B,C;
    
    loop(j,3) {
      p.comp[j] = input[(c-1)*3+j] - input[(b-1)*3+j];
      q.comp[j] = input[(d-1)*3+j] - input[(c-1)*3+j];
      A.comp[j] = input[(a-1)*3+j] - 2*input[(b-1)*3+j] + input[(c-1)*3+j];
      B.comp[j] = input[(c-1)*3+j] - 2*input[(d-1)*3+j] + input[(e-1)*3+j];
    }

    double pp=0,qq=0;
    loop(j,3) {
      pp += p.comp[j]*p.comp[j];
      qq += q.comp[j]*q.comp[j];
    }
    
    double C2 = 0, dotCB = 0, dotCA = 0;
    loop(j,3) {
      C.comp[j] = qq*A.comp[j]-pp*B.comp[j];
      
      dotCB += C.comp[j]*B.comp[j];
      dotCA += C.comp[j]*A.comp[j];
      
      C2 += C.comp[j]*C.comp[j];
    }
    
    result[i] = C2; 
    
    if(grad) {
      loop(j,3) {
        grad[i*n + (a-1)*3+j] =  2*qq*C.comp[j];
        grad[i*n + (e-1)*3+j] =  -2*pp*C.comp[j];
        
        grad[i*n + (b-1)*3+j] = -4*qq*C.comp[j] + 4*p.comp[j]*dotCB;
        grad[i*n + (d-1)*3+j] = 4*pp*C.comp[j] + 4*q.comp[j]*dotCA;
        
        grad[i*n + (c-1)*3+j] =  2*(qq-pp)*C.comp[j] - 4*q.comp[j]*dotCA -4*p.comp[j]*dotCB;
      }
    }
  }
}

void fibnetwork::minimize()
{
  double *x = new double[np*3];
  double *df = new double[np*3];
  loop(i,np) loop(j,3) x[i*3+j] = CPx[i].comp[j];
  
  nlopt_opt opt;

  //establish sizes
  unsigned n = np*3; //number of decision variables
  unsigned m_in = 0; //number of inequality constraints

  //bounds for decision variables
  //double lb[] = { 0.3, -HUGE_VAL, -HUGE_VAL,  -HUGE_VAL, -HUGE_VAL }; /* lower bounds */
  //double ub[] = { HUGE_VAL, HUGE_VAL, 5, HUGE_VAL, HUGE_VAL }; /* lower bounds */
  
  double *lb = new double[np*3];
  double *ub = new double[np*3];
  
  loop(i,np) loop(j,3) {
    lb[i*3+j] = bounds[j][0];
    ub[i*3+j] = bounds[j][1];
  }
  
  //opt = nlopt_create(NLOPT_LN_COBYLA, n);
  //opt = nlopt_create(NLOPT_LD_LBFGS, n);
  opt = nlopt_create(NLOPT_LD_SLSQP, n);
  
  nlopt_set_lower_bounds(opt, lb);
  nlopt_set_upper_bounds(opt, ub);

  nlopt_set_min_objective(opt, VdV, this);
  
  //double tol_eq[]={1e-8,1e-8};
  unsigned m_eq = fixed_dof.size(); //number of equality constraints
  double *fctol = new double[m_eq];  loop(i,m_eq) fctol[i] = 1e-8;
  nlopt_add_equality_mconstraint(opt, m_eq, fixedconstraints, this, fctol);

  unsigned m_eq2 = C1_cps.size(); //number of equality constraints for C1
  double *cctol = new double[m_eq2];  loop(i,m_eq2) cctol[i] = 1e-8;
  nlopt_add_equality_mconstraint(opt, m_eq2, C1_continuity, this, cctol); //C1
  nlopt_add_equality_mconstraint(opt, m_eq2, C2_continuity, this, cctol); //C2
  
  cfac = 0.0;
  nlopt_set_xtol_rel(opt, 1e-6);
  
  double minf;
  
  if(1) {
    nlopt_result status = nlopt_optimize(opt, x, &minf);
    if (status < 0) {
      printf("nlopt failed! Error code: %d\n",status);
      std::cout << "Optimization result: " << nlopt_result_to_string(status) << std::endl;
    }
    else
      printf("found minimum at %0.10g\n", minf);
      loop(i,np) loop(j,3) CPx[i].comp[j] = x[i*3+j];
      printlammps_cps((char*) "Output_cps.lammpstrj",(char*) "a");
      printlammps((char*) "Output.lammpstrj",(char*) "a", NBEZ);  
  }
  
  nlopt_destroy(opt);
  //double fn_energy = VdV(np, x, df, this);
  //printf("FN Energy: %lf\nFN Gradients: \n",fn_energy);
  //loop(i,np) printf("%.4e %.4e %.4e\n",df[i*3],df[i*3+1],df[i*3+2]);
  
  //printf("Fixed dofs: \n");
  //for(std::map<int,double>::iterator it = fixed_dof.begin(); it != fixed_dof.end(); it++)
    //printf("%d %lf\n",it->first,it->second);
    
  //printf("Collinear cps: \n");
  //loop(i,C1_cps.size()) {
    //int a = C1_cps[i].comp[0];
    //int b = C1_cps[i].comp[1];
    //int c = C1_cps[i].comp[2];
    //printf("%d %d %d \n",a,b,c);
  //}
  
  delete [] lb;
  delete [] ub;
  delete [] x;
  delete [] df;
  delete [] fctol;
}

void fibnetwork::assemble_matrices()
{ 
  gsl_matrix_set_zero(GM);
  gsl_matrix_set_zero(GMdot);
  gsl_matrix_set_zero(GN);
  
  double val;
  loop(i,nb) {
    loop(m,4) loop(n,4) {
      
      val = gsl_matrix_get(GM, Bz_list[i].CP[m]-1 , Bz_list[i].CP[n]-1) + gsl_matrix_get(Bz_list[i].M, m,n);
      gsl_matrix_set(GM, Bz_list[i].CP[m]-1 , Bz_list[i].CP[n]-1, val );
      
      val = gsl_matrix_get(GMdot, Bz_list[i].CP[m]-1 , Bz_list[i].CP[n]-1) + gsl_matrix_get(Bz_list[i].Mdot, m,n);
      gsl_matrix_set(GMdot, Bz_list[i].CP[m]-1 , Bz_list[i].CP[n]-1, val );
      
      val = gsl_matrix_get(GN, Bz_list[i].CP[m]-1 , Bz_list[i].CP[n]-1) + gsl_matrix_get(Bz_list[i].N, m,n);
      gsl_matrix_set(GN, Bz_list[i].CP[m]-1 , Bz_list[i].CP[n]-1, val );      
    }
    //gsl_matrix_print(Bz_list[i].M,"Mi");
    //gsl_matrix_print(GM,"GMi");
  }
}

double fibnetwork::dMdQ(int r, int l, int m, int n)
{
  
  double res1 = 0;
  std::vector<int>::iterator it;
  
  for(it = CP_bid[m].begin(); it != CP_bid[m].end(); it++ )    {
    int bid = *it;
    int start = Bz_list[bid-1].CP[0]-1;
    int end = Bz_list[bid-1].CP[3]-1;
    int rank = m-start;
    
    if(r>=start && r<=end && l>=start && l<=end) {
      //printf("\n %d %d %d %d:  %d %d %d %d", r, l, bid, start, r-start, l-start, rank, n);
      //printf("\n %d %d %d %d: (%d) %d %d %d %d", r, l, m, n, bid-1, r-start, l-start, rank, n);
      res1 += Bz_list[bid-1].dM[r-start][l-start][rank][n];
    }
  }

  //return 0.0;
  return res1;
}

void fibnetwork::compute_dH(int flag)
{  
  assemble_matrices(); 
  
  //gsl_matrix_print(GM,"GM0:");
  
  gsl_matrix *GMtmp = gsl_matrix_alloc (np, np);
  gsl_matrix_memcpy(GMtmp, GM);
  
  gsl_permutation *p = gsl_permutation_alloc(np); int s;
  gsl_linalg_LU_decomp(GMtmp, p, &s); gsl_linalg_LU_invert(GMtmp, p, GMinv); 
  
  gsl_permutation_free(p); gsl_matrix_free(GMtmp);

  
  //gsl_matrix_print(GM,"GM1:");
  //gsl_matrix_print(GMinv,"GMinv:");
  //gsl_matrix_print(Wsys,"Wsys:");
  
  gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, GMinv, Wsys, 0.0, pH); //compute dH/dP(mn) = Qdot(mn)
  
  //gsl_matrix_print(pH,"pH:");

  if(1) {
    loop(i,np) loop(j,3) if(!isnan(gsl_matrix_get(pH,i,j))) CPv[i].comp[j] = gsl_matrix_get(pH,i,j);
    
    loop(i,nb) {
      loop(j,4) {
        int cpid = Bz_list[i].CP[j]-1;
        loop(k,3) {
          Bz_list[i].v[j*3+k] = CPv[cpid].comp[k];
        }
      }
      
      //Bz_list[i].momentum();
    }
  }
  
  if(flag) {
    compute_ef();
  
    //gsl_matrix_print(pH,"pH:");
  
    //loop(i,np) loop(j,3) gsl_matrix_set(Fsys, i,j, CPf[i].comp[j]);
    
    //loop(i,np) printf("computeK,CPf: %d %lf %lf %lf\n",i+1,CPf[i].comp[0],CPf[i].comp[1],CPf[i].comp[2]);
    if(1) {
      //gsl_permutation *p = gsl_permutation_alloc(np); int s;
      //gsl_linalg_LU_decomp(GM, p, &s); gsl_linalg_LU_invert(GM, p, GMinv); gsl_permutation_free(p);
      
      //gsl_matrix_print(GMinv,"GMinv:");
      //gsl_matrix_print(Wsys,"Wsys:");
      
      //gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, GMinv, Wsys, 0.0, pH); //compute dH/dP(mn)
      
      //gsl_matrix_print(pH,"pH:");
      
      double kinE1 = 0.0;
      loop(m,np) loop(n,3) {
        kinE1 += 0.5*gsl_matrix_get(pH,m,n)*gsl_matrix_get(Wsys,m,n);
      } 
      kinE = kinE1;
      
      //gsl_matrix_print(pH,"pH:");
      //gsl_matrix_print(Wsys,"Wsys:");
      
      //printf("kinE1: %lf\n",kinE1);
      
      //compute dH/dQ(mn)
      gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, pH, pH, 1.0, GN);
    
      if(1) loop(m,np) loop(n,3) {  
        double dHdQ2=0;
        
        loop(r,np) loop(l,np) {
          dHdQ2 += -0.5*gsl_matrix_get(GN,r,l)*dMdQ(r,l,m,n); //smart mapping required here
        }
        
        gsl_matrix_set(qH, m,n, -CPf[m].comp[n] + dHdQ2);
      }
      
      //gsl_matrix_print(qH,"qH:");
  
    }
  }
}

void fibnetwork::computeK(gsl_matrix *K)
{ 
  assemble_matrices(); 
  
  //loop(i,np) printf("computeK,CPf1: %d %lf %lf %lf\n",i+1,CPf[i].comp[0],CPf[i].comp[1],CPf[i].comp[2]);
  
  
  compute_ef();
  
  loop(i,np) loop(j,3) gsl_matrix_set(Fsys, i,j, CPf[i].comp[j]);
  
  //loop(i,np) printf("computeK,CPf: %d %lf %lf %lf\n",i+1,CPf[i].comp[0],CPf[i].comp[1],CPf[i].comp[2]);
  
  gsl_matrix *GMtmp = gsl_matrix_alloc (np, np);
  gsl_matrix_memcpy(GMtmp, GM);
  
  gsl_permutation *p = gsl_permutation_alloc(np); int s;
  gsl_linalg_LU_decomp(GMtmp, p, &s); gsl_linalg_LU_invert(GMtmp, p, GMinv); 
  
  gsl_permutation_free(p); gsl_matrix_free(GMtmp);
  
  gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, GN, Psys, 1.0, Fsys);
  gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, -1.0, GMdot, Qsys, 1.0, Fsys);
  
  loop(i,np) loop(j,3) gsl_matrix_set(K, i,j, gsl_matrix_get(Qsys, i,j) );
  
  gsl_matrix_scale(Qsys, -RKdamp);	gsl_matrix_add(Fsys, Qsys);
  
  gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, GMinv, Fsys, 0.0, Rsys);
  loop(i,np) loop(j,3) gsl_matrix_set(K, np+i,j, gsl_matrix_get(Rsys, i,j) );

}


/*
void fibnetwork::update_sys(double dh, gsl_matrix *K0, gsl_matrix *K, int c2_flag)
{
  //update Psys, Qsys, CPx, CPv
  loop(i,np) loop(j,3) gsl_matrix_set(Psys, i,j, gsl_matrix_get(K0, i,j) + dh*gsl_matrix_get(K, i,j) );
  loop(i,np) loop(j,3) gsl_matrix_set(Qsys, i,j, gsl_matrix_get(K0, np+i,j) + dh*gsl_matrix_get(K, np+i,j) );
  
  loop(i,np) loop(j,3) CPx[i].comp[j] = gsl_matrix_get(Psys, i,j);
  loop(i,np) loop(j,3) CPv[i].comp[j] = gsl_matrix_get(Qsys, i,j);
  
  if(c2_flag == 1)
    loop(i,C2_cps.size()) {
      int a = C2_cps[i].comp[0];
      int b = C2_cps[i].comp[1];
      int c = C2_cps[i].comp[2];
      int d = C2_cps[i].comp[3];
      int e = C2_cps[i].comp[4];
      
      Vector3 db(CPx[d-1].comp[0]-CPx[b-1].comp[0], CPx[d-1].comp[1]-CPx[b-1].comp[1] , CPx[d-1].comp[2]-CPx[b-1].comp[2] );
      Vector3 eb(CPx[e-1].comp[0]-CPx[b-1].comp[0], CPx[e-1].comp[1]-CPx[b-1].comp[1] , CPx[e-1].comp[2]-CPx[b-1].comp[2] );
      Vector3 ad(CPx[a-1].comp[0]-CPx[d-1].comp[0], CPx[a-1].comp[1]-CPx[d-1].comp[1] , CPx[a-1].comp[2]-CPx[d-1].comp[2] );
      
      double m = sqrt(cross_mag(db,eb)/cross_mag(db,ad));
      double x = 1.0/(1.0+m);
      
      //EDIT to avoid singularities
      if(x<0.1) x = 0.1; if(x>0.9) x = 0.9;
      
      //printf("%d %d %d %d %d %d %lf %lf\n",i,a,b,c,d,e,m,x);
      loop(j,3) CPx[c-1].comp[j] = x*CPx[d-1].comp[j] + (1.0-x)*CPx[b-1].comp[j];
    }
  
  loop(i,nb) {    
    loop(j,4) {
      int cpid = Bz_list[i].CP[j]-1;
      loop(k,3) {
        Bz_list[i].x[j*3+k] = CPx[cpid].comp[k]; Bz_list[i].v[j*3+k] = CPv[cpid].comp[k];
      }
    }
    Bz_list[i].update_bezier();
  }
}

void fibnetwork::integrate_runge_kutta_4(int nsteps, double h)
{
  Psys = gsl_matrix_alloc (np, 3);
  Qsys = gsl_matrix_alloc (np, 3);
  Rsys = gsl_matrix_alloc (np, 3);
  Fsys = gsl_matrix_alloc (np, 3);
  
  gsl_matrix *K0 = gsl_matrix_alloc (2*np, 3);
  gsl_matrix *K1 = gsl_matrix_alloc (2*np, 3);
  gsl_matrix *K2 = gsl_matrix_alloc (2*np, 3);
  gsl_matrix *K3 = gsl_matrix_alloc (2*np, 3);
  gsl_matrix *K4 = gsl_matrix_alloc (2*np, 3);
  gsl_matrix *Kavg = gsl_matrix_alloc (2*np, 3);
  
  gsl_matrix *Minv = gsl_matrix_alloc (np, np);
  
  loop(i,np) loop(j,3) {
    gsl_matrix_set(Psys, i,j, CPx[i].comp[j]);
    gsl_matrix_set(Qsys, i,j, CPv[i].comp[j]);
  }
  
  printf("Time AxialE BendE CohE kinE TotE\n");
  
  for(int step=0; step<nsteps; step++)
  {
    //gsl_matrix_memcpy(Psys0,Psys);
    //gsl_matrix_memcpy(Qsys0,Qsys);
    
    loop(i,np) loop(j,3) gsl_matrix_set(K0, i,j, gsl_matrix_get(Psys, i,j));
    loop(i,np) loop(j,3) gsl_matrix_set(K0, np+i,j, gsl_matrix_get(Qsys, i,j));
    
    computeK(K1); update_sys(0.5*h,K0,K1,0); //RK step 1
    computeK(K2); update_sys(0.5*h,K0,K2,0); //RK step 2
    computeK(K3); update_sys(h,K0,K3,0); //RK step 3
    
    computeK(K4); //RK step 4
    loop(i,2*np) loop(j,3) {
      double Kij = gsl_matrix_get(K1, i,j) + 2*gsl_matrix_get(K2, i,j) + 2*gsl_matrix_get(K3, i,j) + gsl_matrix_get(K4, i,j) ;
      gsl_matrix_set(Kavg, i,j, Kij);
    }
    update_sys((h/6.0),K0,Kavg,1); 
    
    {
      assemble_matrices();
      compute_ef();
      printlammps_cps((char*) "Output_cps.lammpstrj",(char*) "a");
      printlammps((char*) "Output.lammpstrj",(char*) "a", NBEZ);
      time = step*h; tstep = step;
      printf("%lf %lf %lf %lf %lf %lf\n",time, axialE, bendE, cohE, kinE, totE);
    }
  }
  
  gsl_matrix_free(Psys);
  gsl_matrix_free(Qsys);
  gsl_matrix_free(Rsys);
  gsl_matrix_free(Fsys);
  
  gsl_matrix_free(K0);  gsl_matrix_free(K1);  gsl_matrix_free(K2);
  gsl_matrix_free(K3);  gsl_matrix_free(K4);  gsl_matrix_free(Kavg);
}
*/
