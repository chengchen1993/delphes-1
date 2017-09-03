/*
 *  Delphes: a framework for fast simulation of a generic collider experiment
 *  Copyright (C) 2012-2014  Universite catholique de Louvain (UCL), Belgium
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */


/** \class MomentumSmearing
 *
 *  Performs transverse momentum resolution smearing.
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *
 */

#include "modules/MomentumSmearing.h"

#include "classes/DelphesClasses.h"
#include "classes/DelphesFactory.h"
#include "classes/DelphesFormula.h"

#include "ExRootAnalysis/ExRootResult.h"
#include "ExRootAnalysis/ExRootFilter.h"
#include "ExRootAnalysis/ExRootClassifier.h"

#include "TMath.h"
#include "TString.h"
#include "TFormula.h"
#include "TRandom3.h"
#include "TObjArray.h"
#include "TDatabasePDG.h"
#include "TLorentzVector.h"

#include <algorithm> 
#include <stdexcept>
#include <iostream>
#include <sstream>

using namespace std;

//------------------------------------------------------------------------------

MomentumSmearing::MomentumSmearing() :
	fFormula(0), fItInputArray(0)
{
	fFormula = new DelphesFormula;
}

//------------------------------------------------------------------------------

MomentumSmearing::~MomentumSmearing()
{
	if(fFormula) delete fFormula;
}

//------------------------------------------------------------------------------

void MomentumSmearing::Init()
{
	// read resolution formula

	fFormula->Compile(GetString("ResolutionFormula", "0.0"));
	ff = GetDouble("fff", 0);
	// import input array

	fInputArray = ImportArray(GetString("InputArray", "ParticlePropagator/stableParticles"));
	fItInputArray = fInputArray->MakeIterator();

	// create output array

	fOutputArray = ExportArray(GetString("OutputArray", "stableParticles"));
}

//------------------------------------------------------------------------------

void MomentumSmearing::Finish()
{
	if(fItInputArray) delete fItInputArray;
}

//------------------------------------------------------------------------------

void MomentumSmearing::Process()
{
	Candidate *candidate, *mother;
	Double_t pt, eta, phi, e, res;

	fItInputArray->Reset();
	while((candidate = static_cast<Candidate*>(fItInputArray->Next())))
	{
		const TLorentzVector &candidatePosition = candidate->Position;
		const TLorentzVector &candidateMomentum = candidate->Momentum;
		eta = candidatePosition.Eta();
		phi = candidatePosition.Phi();
		pt = candidateMomentum.Pt();
		e = candidateMomentum.E();
		//cout<<ff<<endl;
		if(ff<1){
			res = fFormula->Eval(pt, eta, phi, e);
		}
		else{

			//TVector3 currP = a_MCP->getMomentum();
			//int currPID = a_MCP->getPDG();
			//float PAmplitude = currP.Mag();
			//float cosTheta = fabs(a_MCP->getMomentum()[2]/PAmplitude);
			float ObjkCoeff = 0; 
			float kCoeffLow = 0;
			float kCoeffHigh = 0; 
			float EnDisToLow = 0;
			float EnDisToHigh = 0; 	
			float Scalefactor = 0; 	//Correct EndCap divergence
			double cosTheta=abs((1-exp(-2*eta))/(1+exp(-2*eta)));

			float RefEnergy[7] = {5, 10, 20, 40, 60, 80, 100};
			float kCoeff[7] = {0.000211212, 0.000107858, 6.01923e-05, 3.69668e-05, 2.8476e-05, 2.44015e-05, 2.24783e-05};	

			if(e < 5)
			{ 
				ObjkCoeff = 0.000267917;
			}
			else if(e > 100)
			{ 
				ObjkCoeff = 3.61945e-05;
			}
			else
			{
				for(int i = 0; i < 6; i++)
				{
					if(e < RefEnergy[i+1] && e > RefEnergy[i])
					{
						kCoeffLow = kCoeff[i];
						kCoeffHigh = kCoeff[i+1];
						EnDisToLow = e - RefEnergy[i];	
						EnDisToHigh = RefEnergy[i+1] - e;
						ObjkCoeff = (kCoeffLow*EnDisToHigh + kCoeffHigh*EnDisToLow)/(RefEnergy[i+1] - RefEnergy[i]);
					}
				}
			}

			if(cosTheta > 0.86)	//Scale as effective R^2
			{
				Scalefactor = (1.0/(cosTheta*cosTheta) - 1)*2.96; //2.96 = (Half_Z/Radius)**2
				ObjkCoeff = ObjkCoeff*1.0/Scalefactor; 
			}

			res = ObjkCoeff*e; 	// To be extended
		}

		res = ( res > 1.0 ) ? 1.0 : res; 
		//Reso = ( Reso > 1.0 ) ? 1.0 : Reso;

		pt = LogNormal(pt, res * pt );

		//   pt = LogNormal(pt, Reso * pt );

		//if(pt <= 0.0) continue;

		mother = candidate;
		candidate = static_cast<Candidate*>(candidate->Clone());
		eta = candidateMomentum.Eta();
		phi = candidateMomentum.Phi();
		candidate->Momentum.SetPtEtaPhiE(pt, eta, phi, pt*TMath::CosH(eta));
		//candidate->TrackResolution = fFormula->Eval(pt, eta, phi, e);
		candidate->TrackResolution = res;
		//candidate->TrackResolution = Reso;
		candidate->AddCandidate(mother);

		fOutputArray->Add(candidate);
	}
}
//----------------------------------------------------------------

Double_t MomentumSmearing::LogNormal(Double_t mean, Double_t sigma)
{
	Double_t a, b;

	if(mean > 0.0)
	{
		b = TMath::Sqrt(TMath::Log((1.0 + (sigma*sigma)/(mean*mean))));
		a = TMath::Log(mean) - 0.5*b*b;

		return TMath::Exp(a + b*gRandom->Gaus(0.0, 1.0));
	}
	else
	{
		return 0.0;
	}
}


//------------------------------------------------------------------------------
