/***************************************************************************
 *
 * Authors:    Carlos Oscar Sanchez Sorzano coss@cnb.csic.es
 *
 * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
 * 02111-1307  USA
 *
 *  All comments concerning this program package may be sent to the
 *  e-mail address 'xmipp@cnb.csic.es'
 ***************************************************************************/

#include "classify_significant.h"
#include <data/mask.h>

// Empty constructor =======================================================
ProgClassifySignificant::~ProgClassifySignificant()
{
	for (size_t i=0; i<projector.size(); i++)
		delete projector[i];
	for (size_t i=0; i<subsetProjections.size(); i++)
		delete subsetProjections[i];
}

// Read arguments ==========================================================
void ProgClassifySignificant::readParams()
{
    fnVols = getParam("--ref");
    fnIds = getParam("--id");
    fnAngles = getParam("--angles");
    fnOut = getParam("-o");
    pad = getIntParam("--padding");
    wmin = getDoubleParam("--minWeight");
}

// Show ====================================================================
void ProgClassifySignificant::show()
{
    if (!verbose)
        return;
    std::cout
    << "Reference volumes:   " << fnVols            << std::endl
    << "Ids:                 " << fnIds             << std::endl
    << "Angles:              " << fnAngles          << std::endl
    << "Output:              " << fnOut             << std::endl
    << "Padding factor:      " << pad               << std::endl
	<< "Min. Weight:         " << wmin              << std::endl
    ;
}

// usage ===================================================================
void ProgClassifySignificant::defineParams()
{
    addUsageLine("Classify a set of images into different classes. See protocol_reconstruct_heterogeneous");
    addParamsLine("   --ref <metadata>            : Reference volumes");
    addParamsLine("   --id <metadata>             : List of itemIds to classified. Sorted.");
    addParamsLine("   --angles <metadata>         : Angles assigned. Each image should have one or several angles");
    addParamsLine("                               : for each volume. The assignment per volume should be organized");
    addParamsLine("                               : in different blocks");
    addParamsLine("   -o <metadata>               : Output metadata with a set of angles per volume");
    addParamsLine("  [--padding <p=2>]            : Padding factor");
    addParamsLine("  [--minWeight <p=0.1>]        : Minimum weight");
}

// Produce side information ================================================
void ProgClassifySignificant::produceSideInfo()
{
	if (verbose>0)
		std::cerr << "Producing side info ..." << std::endl;
    // Read the reference volumes
    Image<double> V;
    MetaData mdVols, mdAngles;
    mdVols.read(fnVols);
    FileName fnVol;
    int i=1;
    FOR_ALL_OBJECTS_IN_METADATA(mdVols)
    {
    	mdVols.getValue(MDL_IMAGE,fnVol,__iter.objId);
    	std::cout << fnVol << std::endl;
        V.read(fnVol);
        V().setXmippOrigin();
        projector.push_back(new FourierProjector(V(),pad,0.5,BSPLINE3));
        currentRowIdx.push_back(0);

        mdAngles.read(formatString("angles_%02d@%s",i,fnAngles.c_str()));
        VMetaData *vmd=new VMetaData();
        mdAngles.asVMetaData(*vmd);
        setAngles.push_back(*vmd);
        classifiedAngles.push_back(*(new VMetaData()));

        subsetAngles.push_back(*(new VMetaData()));
        subsetProjectionIdx.push_back(* (new std::vector<size_t>));
        i+=1;
    }

    // Read the Ids
    MetaData mdIds;
    mdIds.read(fnIds);
    mdIds.getColumnValues(MDL_PARTICLE_ID,setIds);
}

//#define DEBUG
void ProgClassifySignificant::generateProjection(size_t volumeIdx, size_t poolIdx, MDRow &currentRow)
{
	double rot, tilt, psi, x, y;
	bool flip;
	Matrix2D<double> A;
	currentRow.getValue(MDL_ANGLE_ROT,rot);
	currentRow.getValue(MDL_ANGLE_TILT,tilt);
	currentRow.getValue(MDL_ANGLE_PSI,psi);
	currentRow.getValue(MDL_SHIFT_X,x);
	currentRow.getValue(MDL_SHIFT_Y,y);
	currentRow.getValue(MDL_FLIP,flip);
	A.initIdentity(3);
	MAT_ELEM(A,0,2)=x;
	MAT_ELEM(A,1,2)=y;
	if (flip)
	{
		MAT_ELEM(A,0,0)*=-1;
		MAT_ELEM(A,0,1)*=-1;
		MAT_ELEM(A,0,2)*=-1;
	}
	projectVolume(*(projector[volumeIdx]), Paux, (int)XSIZE(Iexp()), (int)XSIZE(Iexp()),  rot, tilt, psi);

	if (poolIdx>=subsetProjections.size())
		subsetProjections.push_back(new MultidimArray<double>);
	subsetProjections[poolIdx]->resizeNoCopy((int)XSIZE(Iexp()), (int)XSIZE(Iexp()));
	applyGeometry(LINEAR,*(subsetProjections[poolIdx]),Paux(),A,IS_INV,DONT_WRAP,0.);

#ifdef DEBUG
	std::cout << "Row: " << " rot: " << rot << " tilt: " << tilt
			  << " psi: " << psi << " sx: " << x << " sy: " << y
			  << " flip: " << flip << std::endl;
	Image<double> save;
	save()=*(subsetProjections[poolIdx]);
	save.write(formatString("PPPtheo%02d.xmp",poolIdx));
#endif
}

// Select the subset associated to a particleId
void ProgClassifySignificant::selectSubset(size_t particleId)
{
	size_t poolIdx=0;
	FileName fnImg;
	for (size_t i=0; i<projector.size(); i++)
	{
		subsetAngles[i].clear();
		subsetProjectionIdx[i].clear();
		size_t crIdx=currentRowIdx[i];
		if (crIdx>=setAngles[i].size())
			return;
		MDRow & currentRow=setAngles[i][crIdx];
		if (i==0) // First time we see this image
		{
			currentRow.getValue(MDL_IMAGE,fnImg);
			Iexp.read(fnImg);
		}
		size_t currentParticleId;
		currentRow.getValue(MDL_PARTICLE_ID,currentParticleId);
		size_t idxMax=setAngles[i].size();
		while (currentParticleId<=particleId)
		{
			if (currentParticleId==particleId)
			{
				subsetAngles[i].push_back(currentRow);
				subsetProjectionIdx[i].push_back(poolIdx);
				generateProjection(i,poolIdx,currentRow);
				poolIdx+=1;
			}
			crIdx+=1;
			if (crIdx<idxMax)
			{
				currentRow=setAngles[i][crIdx];
				currentRow.getValue(MDL_PARTICLE_ID,currentParticleId);
			}
			else
				break;
		}
		currentRowIdx[i]=crIdx;
	}
#ifdef DEBUG
	std::cout << "Reading " << fnImg << std::endl;
	char c;
	std::cout << "Press any key" << std::endl;
	std::cin >> c;
#endif
}
#undef DEBUG

void computeWeightedCorrelation(const MultidimArray<double> &I1, const MultidimArray<double> &I2, const MultidimArray<double> &Iexp,
		const MultidimArray<double> &Idiff, double &corr1exp, double &corr2exp)
{
	double mean, std;
	Idiff.computeAvgStdev(mean,std);
	double threshold=std;

	corr1exp=corr2exp=0.0;

	// Estimate the mean and stddev within the mask
	double N=0;
	double sumWI1=0, sumWI2=0, sumWIexp=0;
	double sumI1=0, sumI2=0, sumIexp=0;
	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Idiff)
	{
		double p1=DIRECT_MULTIDIM_ELEM(I1,n);
		double p2=DIRECT_MULTIDIM_ELEM(I2,n);
		double pexp=DIRECT_MULTIDIM_ELEM(Iexp,n);
		sumI1+=p1;
		sumI2+=p2;
		sumIexp+=pexp;
		if (DIRECT_MULTIDIM_ELEM(Idiff,n)>threshold)
		{
			sumWI1+=p1;
			sumWI2+=p2;
			sumWIexp+=pexp;
			N+=1.0;
		}
	}

	// Estimate the weighted correlation
	if (N>0)
	{
		double iN=1.0/N;
		double isize=1.0/MULTIDIM_SIZE(Idiff);
		double avgW1=sumWI1*iN;
		double avgW2=sumWI2*iN;
		double avgWExp=sumWIexp*iN;

		double sumWI1exp=0.0, sumWI2exp=0.0, sumWI1I1=0.0, sumWI2I2=0.0, sumWIexpIexp=0.0;
		double sumI1exp=0.0,  sumI2exp=0.0,  sumI1I1=0.0,  sumI2I2=0.0,  sumIexpIexp=0.0;
		FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Idiff)
		{
			double p1=DIRECT_MULTIDIM_ELEM(I1,n);
			double p2=DIRECT_MULTIDIM_ELEM(I2,n);
			double pexp=DIRECT_MULTIDIM_ELEM(Iexp,n);
			double p1a=p1-avgW1;
			double p2a=p2-avgW2;
			double pexpa=pexp-avgWExp;
			sumI1exp+=p1a*pexpa;
			sumI2exp+=p2a*pexpa;
			sumI1I1 +=p1a*p1a;
			sumI2I2 +=p2a*p2a;
			sumIexpIexp +=pexp*pexp;

			if (DIRECT_MULTIDIM_ELEM(Idiff,n)>threshold)
			{
				p1a=p1-avgW1;
				p2a=p2-avgW2;
				pexpa=pexp-avgWExp;

				double w=DIRECT_MULTIDIM_ELEM(Idiff,n);
				double wp1a=w*p1a;
				double wp2a=w*p2a;
				double wpexpa=w*pexpa;

				sumWI1exp+=wp1a*pexpa;
				sumWI2exp+=wp2a*pexpa;
				sumWI1I1 +=wp1a*p1a;
				sumWI2I2 +=wp2a*p2a;
				sumWIexpIexp +=wpexpa*pexpa;
			}
		}

		sumWI1exp*=iN;
		sumWI2exp*=iN;
		sumWI1I1*=iN;
		sumWI2I2*=iN;
		sumWIexpIexp*=iN;

		sumI1exp*=isize;
		sumI2exp*=isize;
		sumI1I1*=isize;
		sumI2I2*=isize;
		sumIexpIexp*=isize;

		double corrW1exp=sumWI1exp/sqrt(sumWI1I1*sumWIexpIexp);
		double corrW2exp=sumWI2exp/sqrt(sumWI2I2*sumWIexpIexp);
		double corrN1exp=sumI1exp/sqrt(sumI1I1*sumIexpIexp);
		double corrN2exp=sumI2exp/sqrt(sumI2I2*sumIexpIexp);
		if (corrW1exp>0 && corrN1exp>0)
			corr1exp=sqrt(corrW1exp*corrN1exp);
		else
			corr1exp=-1;
		if (corrW2exp>0 && corrN2exp>0)
			corr2exp=sqrt(corrW2exp*corrN2exp);
		else
			corr2exp=-1;
	}
}

void ProgClassifySignificant::updateClass(int n, double wn)
{
	double CCbest=-1e38;
	int iCCbest=-1;
	VMetaData &subsetAngles_n=subsetAngles[n];
	for (int i=0; i<subsetAngles_n.size(); i++)
	{
		double cc;
		subsetAngles_n[i].getValue(MDL_MAXCC,cc);
		if (cc>CCbest)
		{
			CCbest=cc;
			iCCbest=i;
		}
	}
	if (iCCbest>=0)
	{
		MDRow newRow=subsetAngles_n[iCCbest];
		// COSS newRow.setValue(MDL_WEIGHT,wn);
		classifiedAngles[n].push_back(newRow);
	}
}

//#define DEBUG
void ProgClassifySignificant::run()
{
	show();
	produceSideInfo();
	MultidimArray<double> Idiff;

	if (verbose>0)
	{
		std::cerr << "Classifying images ..." << std::endl;
		init_progress_bar(setIds.size());
	}

	Matrix1D<double> winning(projector.size());
	Matrix1D<double> corrDiff(projector.size());
	Matrix1D<double> weight;
	for (size_t iid=0; iid<setIds.size(); iid++)
	{
		size_t id=setIds[iid];
		selectSubset(id);
		const MultidimArray<double> &mIexp=Iexp();
		winning.initZeros();
		corrDiff.initZeros();
		for (size_t ivol1=0; ivol1<projector.size(); ivol1++)
		{
			std::vector<size_t> &subset1=subsetProjectionIdx[ivol1];
			for (size_t ivol2=ivol1+1; ivol2<projector.size(); ivol2++)
			{
				std::vector<size_t> &subset2=subsetProjectionIdx[ivol2];
				for (size_t i1=0; i1<subset1.size(); i1++)
				{
					MultidimArray<double> &I1=*(subsetProjections[subset1[i1]]);
					for (size_t i2=0; i2<subset2.size(); i2++)
					{
						MultidimArray<double> &I2=*(subsetProjections[subset2[i2]]);
						Idiff=I1;
						Idiff-=I2;
						Idiff.selfABS();

						double corr1exp, corr2exp;
						computeWeightedCorrelation(I1, I2, mIexp, Idiff, corr1exp, corr2exp);
						double corrDiff12=corr1exp-corr2exp;
						if (corrDiff12>0 && corr1exp>0)
						{
							VEC_ELEM(winning,ivol1)+=1;
							VEC_ELEM(corrDiff,ivol1)+=corrDiff12;
							VEC_ELEM(corrDiff,ivol2)-=corrDiff12;
						}
						else if (corrDiff12<0 && corr2exp>0)
						{
							VEC_ELEM(winning,ivol2)+=1;
							VEC_ELEM(corrDiff,ivol2)-=corrDiff12;
							VEC_ELEM(corrDiff,ivol1)+=corrDiff12;
						}

//						Image<double> save;
//						save()=Iexp();
//						save.write("PPPexp.xmp");
//						save()=I1;
//						save.write("PPP1.xmp");
//						save()=I2;
//						save.write("PPP2.xmp");
//						save()=Idiff;
//						save.write("PPPdiff.xmp");
//						std::cout << "corr1exp=" << corr1exp << " corr2exp=" << corr2exp << std::endl;
//						char c;
//						std::cout << "Press any key" << std::endl;
//						std::cin >> c;
					}
				}
			}
		}
		double iNcomparisons=1.0/winning.sum();
		winning*=iNcomparisons;
		weight=corrDiff;
		weight*=iNcomparisons;
		weight*=winning;

//		std::cout << corrDiff << std::endl;
//		std::cout << winning << std::endl;
//		std::cout << weight << std::endl;
//		char c;
//		std::cout << "Press any key" << std::endl;
//		std::cin >> c;

		int nBest=weight.maxIndex();
		double wBest=VEC_ELEM(weight,nBest);
		if (wBest>0)
			updateClass(nBest,wBest);
		if (verbose>0)
			progress_bar(iid);
	}
	progress_bar(setIds.size());

	// Write output
	MetaData md;
	for (size_t ivol=0; ivol<projector.size(); ivol++)
	{
		size_t objId=md.addObject();
		md.setValue(MDL_REF3D,(int)ivol+1,objId);
		md.setValue(MDL_CLASS_COUNT,classifiedAngles[ivol].size(),objId);
	}
	md.write("classes@"+fnOut);
	for (size_t ivol=0; ivol<projector.size(); ivol++)
	{
		md.clear();
		md.fromVMetaData(classifiedAngles[ivol]);
		double currentWmax=md.getColumnMax(MDL_WEIGHT);
		double currentWmin=md.getColumnMin(MDL_WEIGHT);
		if (currentWmax>currentWmin)
			md.operate(formatString("weight=%f*(weight-%f)+%f",(1.0-wmin)/(currentWmax-currentWmin),currentWmin,wmin));
		else
			md.operate(formatString("weight=%f",wmin));
		md.setValueCol(MDL_REF3D,(int)ivol+1);
		md.write(formatString("class%06d_images@%s",ivol+1,fnOut.c_str()),MD_APPEND);
	}
}
#undef DEBUG
