/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   See http://www.plumed-code.org for more information.

   plumed is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   plumed is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with plumed.  If not, see <http://www.gnu.org/licenses/>.

   CmuMD method file, 
   see (and cite) Perego, Salvalaglio, Parrinello J. Chem. Phys. 142 144113 (2015) 
   http://scitation.aip.org/content/aip/journal/jcp/142/14/10.1063/1.4917200

   CmuMD, cannibalistic (growth & dissolution)
   see (and cite) Karmakar, Piaggi, Perego, Parrinello J. Chem. Theory Comput., 2018, 14 (5), pp 2678–2683 
   https://pubs.acs.org/doi/abs/10.1021/acs.jctc.8b00191

   CmuMD, nucleation
   see (and cite) Karmakar, Piaggi, Parrinello 
   Molecular Dynamics Simulations of Crystal Nucleation from Solution at Constant Chemical Potential
   https://arxiv.org/abs/1907.04037
   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
 
#include "colvar/Colvar.h"
#include "colvar/ActionRegister.h"

#include <string>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>

using namespace std;

namespace PLMD{
  namespace colvar{

    //+PLUMEDOC COLVAR NSHELLSP 
    /*
      Calculates the solute and solvent concentration in a spherical SHELL defined by inner and outer radius
      The solute/solvent molecules are restrained in that SHELL, Contron Region (CR)
      example: 
           NSHELLSP ...
           LABEL=n_na
           GROUP=na
           NASV=1
           NST=1000
           NAST=1
           DIN=1.7
           DOUT=2.5
           DF=2.6
           COIN=1.7 COOUT=2.5
           WF=0.2 WOUT=0.05 WIN=0.05
           NOSCALE
         ... NSHELLSP

         RESTRAINT ARG=n_na AT=5.0 KAPPA=10000.0 LABEL=nares
    */
    //+ENDPLUMEDOC
   
    class NshellSp : public Colvar {
      bool issolute;
      bool isnotscaled;
      bool isFirstStep;
      int storeHalfBin;
      int  Na_sv_permol, Na_st_permol, N_sv, N_st, N_mol;
      double  D_IN, D_OUT, D_CR, D_F;
      double  w_force, w_in, w_out, co_out, co_in, co_f;
      
    public:
      NshellSp(const ActionOptions&);
      virtual void calculate();
      static void registerKeywords(Keywords& keys);
      double sigmon(double d, double Coff);
      double sigmoff(double d, double Coff);
      double dsig(double d, double Coff);
      ofstream fdbg;
    };

    PLUMED_REGISTER_ACTION(NshellSp,"NSHELLSP")

    void NshellSp::registerKeywords(Keywords& keys){
      Colvar::registerKeywords(keys);
      keys.add("atoms","GROUP","the group of atoms involved in the calculation");
      keys.add("optional","NAST","Number of atoms per Solute molecule");
      keys.add("optional","NST","Solute tot atoms");
      keys.add("optional","NASV","Number of atoms per Solvent molecule");
      keys.add("optional","NSV","Solvent tot atoms");
      keys.add("compulsory","DIN","inner CR boundary, bias region (BR)");
      keys.add("compulsory","DOUT","outer CR boundary");
      keys.add("optional","DF","Force radius");
      keys.add("compulsory","WF","force sigma length");
      keys.add("optional","COF","force sigma cutoff");
      keys.add("optional","WIN","in sigma length");
      keys.add("optional","COIN","in sigma cutoff");
      keys.add("optional","WOUT","out sigma length");
      keys.add("optional","COOUT","out sigma cutoff");

      keys.addFlag("NOSCALE",false,"use absolute length units");
      keys.remove("NOPBC");
    }

    NshellSp::NshellSp(const ActionOptions&ao):
      PLUMED_COLVAR_INIT(ao),
      //init bool parameters
      isnotscaled(false)
    {
      
      //Read atom group
      vector<AtomNumber> at_list;
      parseAtomList("GROUP",at_list);

      // Solvent atoms
      Na_sv_permol=1; //default
      parse("NASV",Na_sv_permol); //number of atoms per Solvent molecule, COM atoms could be used for molecules
      parse("NSV",N_sv); //get total number of Solvent atoms

      // Solute atoms
      N_st=0; //default
      Na_st_permol=1;
      parse("NAST",Na_st_permol); //number of atoms per Solute molecule, COM atoms could be used externally
      parse("NST",N_st); //get total numner of Solute atoms
      
      //Solution atoms
      N_sv=at_list.size()-N_st; //Number of solvent atoms, one atom description
      N_mol=N_sv+N_st; //Number of total molecules
      
      //log.printf("Number of atoms:\tw %d\tu %d\n",N_sv,N_st);
      log.printf("Number of atoms:\ttot %d\t w %d\tu %d\n",N_mol,N_sv,N_st);

      //Parameters (force position and switching function temperature)
      parse("DIN",D_IN);      //inner CR boundary
      parse("DOUT",D_OUT);    //outer CR boundary

      D_F=D_OUT;        //initialize D_F, force just outside of the CR boundary                        
      D_CR=D_OUT-D_IN;

      parse("DF",D_F); // user-defined force radius 
      if(D_F<D_OUT){   //re-initialize D_F if inside CR
	D_F=D_OUT;     // place it at the outer boundary
	log.printf("D_F inside CR region, reset at the boundary");
      }
   
     // Parameters for the switching functions
      parse("WF",w_force); //Fermi Fun T at DF      
      co_f=2.0; //initialize cut-off in
      parse("COF",co_f); //cut-off for Fermi f

      w_in=w_force; //initialize w_in
      parse("WIN",w_in); //Fermi Fun T at CRin
      co_in=co_f; //initialize cut-off in
      parse("COIN",co_in); //cut-off for Fermi f

      w_out=w_force; //initialize w_out
      parse("WOUT",w_out); //Fermi Fun T at CRout
      co_out=co_f; //initialize cut-off in
      parse("COOUT",co_out); //cut-off for Fermi f
      
      log.printf("Geometry:\tD_IN %lf\tD_OUT %lf\tD_CR %lf\tD_F %lf\n",D_IN,D_OUT,D_CR,D_F);
      log.flush();
      
      //other bool parameters 
      parseFlag("NOSCALE",isnotscaled);
      checkRead();
      addValueWithDerivatives(); 
      setNotPeriodic();
      
      //log atom lists
      log.printf("  of %d atoms\n",at_list.size());
      for(unsigned int i=0;i<at_list.size();++i){
	log.printf("  %d", at_list[i].serial());
      }
      log.printf("  \n");
      if(N_st>0){ 
	log.printf("of which the first %d are solute atoms\n",N_st);
      }
      requestAtoms(at_list);
      log.printf("  \n");
      log.flush();       
      isFirstStep=true;
    }

// Switching functions
  double NshellSp::sigmon(double d, double Coff){
    double sig;
    if( d < Coff){
      sig=0.0;
    }else{
      sig=1.0/(exp(-d)+1.0);
    }
    return(sig);
  }


  double NshellSp::sigmoff(double d, double Coff){
    double sig;
    if(d > Coff){
      sig=0.0;
    }else{
      sig=1.0/(exp(d)+1.0);
    }
    return(sig);
  }
 

 // Bell-shaped function (G) in the Force equation 
  double NshellSp::dsig(double d, double Coff){
    double dsig;
    if(fabs(d) > Coff){
      dsig=0.0;
    }else{
      dsig=0.5/(1.0+cosh(d));
    }
    return(dsig);
  }

    
 // calculator
  void NshellSp::calculate()    
  {
      double n_CR;
      Tensor virial;
      virial.zero();  //no virial contribution, required for NPT simulation, will be added in the next version
                      // for NVT simulation it is not necessary (restrained solution conc. is enough)
      
    //Vector deriv;
      vector<Vector> deriv(getNumberOfAtoms());

      Vector ze;
      ze.zero();
      fill(deriv.begin(), deriv.end(), ze); //initialize derivatives

    //Parallel parameters
      unsigned int stride;
      unsigned int rank; 
      stride=comm.Get_size();  //Number of processes
      rank=comm.Get_rank();    //Rank of present process
      
    //Box size
      double LBC[3];
      for(int i=0;i<3;++i) LBC[i]=getBox()[i][i]; 


    //CR volume 
      double VCR;
      VCR=(4.0/3.0)*pi*(D_OUT*D_OUT*D_OUT-D_IN*D_IN*D_IN); //CR volume 
    //log.printf("check distance: \tDIN %lf\tDOUT, %lf\tDCR, %lf\tVCR %lf\n",D_IN,D_OUT,D_CR,VCR);
      
  //Evaluate concentration and derivatives
    //if isolute is true C counts the solute molecules, else the solvent ones
      n_CR=0.0;
      double d,din,dout,dF,n_rx,dfunc; 

    //Solute and Solvent position matrix allocation 
      vector<Vector> solv_x(N_sv);        // [i] x y z
      vector<Vector> solut_x(N_st);

      int k;
      if(N_st == 0){  //if solvent species is restrained
        for(int i=rank; i<N_sv; i+=stride){
           dfunc=0;
           n_rx=0;
           vector<Vector> id(N_sv);
         
           solv_x[i].zero();
           solv_x[i]=getPosition(N_st+i*Na_sv_permol);
           
        // Distance of each molecule from the box center
           id[i].zero();
           id[i] = pbcDistance(solv_x[i],Vector(LBC[0]/2,LBC[1]/2,LBC[2]/2)); 
           d = id[i].modulo();
           double inv_d=1.0/d;

           din=(d-D_IN)/w_in;
           dout=(d-D_OUT)/w_out;
        // sigmaon at din, sigmaoff at dout
           n_rx=sigmon(din,co_in)*sigmoff(dout,co_out); 

           dF=(d-D_F)/w_force;
           dfunc=(dsig(dF,co_f))/w_force;   

           n_CR+=n_rx; //update CV (for now this is the number of molcules)
          
         //include the derivatives here
           k=N_st+i*Na_sv_permol; //atom counter 
           for (unsigned int ix=0; ix<3; ix++) {          
             deriv[k][ix]  +=  dfunc*inv_d*id[k][ix]/VCR; // Derivative of dzin/dx, dzin/dy, & dzin/dz
           }
         }
	vector<Vector>().swap(solv_x);
      }else{   //if solute species is restrained
        for(int i=rank; i<N_st; i+=stride){
           dfunc=0;
           n_rx=0;
           vector<Vector> id(N_st);
         
           solut_x[i].zero();
           solut_x[i]=getPosition(i*Na_st_permol);
        
        // Distance of each molecule from the box center
           id[i].zero();
           id[i] = pbcDistance(solut_x[i],Vector(LBC[0]/2,LBC[1]/2,LBC[2]/2)); 
           d = id[i].modulo();
           double inv_d=1.0/d;

           din=(d-D_IN)/w_in;
           dout=(d-D_OUT)/w_out;
           n_rx=sigmon(din,co_in)*sigmoff(dout,co_out); 

           dF=(d-D_F)/w_force;
           dfunc=(dsig(dF,co_f))/w_force; 
         
           n_CR+=n_rx; //update CV (for now this is the number of molcules)
         
         //include the derivatives here
           k=i*Na_st_permol; 
           for (unsigned int ix=0; ix<3; ix++) {         
             deriv[k][ix]  +=  dfunc*inv_d*id[k][ix]/VCR; // Derivative of dzin/dx, dzin/dy, & dzin/dz
           }
	}
	vector<Vector>().swap(solut_x);
      }
      
      comm.Sum(deriv);
      comm.Sum(n_CR);
      comm.Sum(virial);

      int Natot=N_st+N_sv;
      for(int i=0; i< Natot; ++i){
        setAtomsDerivatives(i, deriv[i]);
      }
      
      setValue(n_CR/VCR);
      setBoxDerivatives(virial);
    //setBoxDerivativesNoPbc();

    }
  }  
}    
