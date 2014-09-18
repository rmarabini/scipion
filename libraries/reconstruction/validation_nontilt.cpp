/***************************************************************************
 * Authors:     AUTHOR_NAME (jvargas@cnb.csic.es)
 *
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

#include "validation_nontilt.h"
#include <cstdlib>      // std::rand, std::srand
#include <ctime>        // std::time
#include <algorithm>
#include "data/sampling.h"


void ProgValidationNonTilt::readParams()
{

    fnIn = getParam("-i");
    fnDir = getParam("--odir");
    fnSym = getParam("--sym");
    fnInit = getParam("--volume");
    alpha0 = getDoubleParam("--alpha0");
    //alphaF = getDoubleParam("--alphaF");
    //Niter = getIntParam("--iter");
    //keepIntermediateVolumes = checkParam("--keepIntermediateVolumes");
    angularSampling=getDoubleParam("--angularSampling");
    sampling_rate = getDoubleParam("--sampling_rate");
    //maxShift=getDoubleParam("--maxShift");
    //tilt0=getDoubleParam("--minTilt");
    //tiltF=getDoubleParam("--maxTilt");
    //useImed=checkParam("--useImed");
    //strict=checkParam("--strictDirection");
    //angDistance=getDoubleParam("--angDistance");
    //Nvolumes=getIntParam("--numberOfVolumes");
    Nvolumes = 1;

}

void ProgValidationNonTilt::defineParams()
{
    //usage
    addUsageLine("Validate a 3D reconstruction from its projections attending to directionality and spread of the angular assignments from a given significant value");
    //params
    addParamsLine("   -i <md_file>                : Metadata file with input projections");
    //addParamsLine("  [--numberOfVolumes <N=1>]    : Number of volumes to reconstruct");
    addParamsLine("  [--volume <md_file=\"\">]    : Volume to validate");
    addParamsLine("  [--odir <outputDir=\".\">]   : Output directory");
    addParamsLine("  [--sym <symfile=c1>]         : Enforce symmetry in projections");
    //addParamsLine("  [--iter <N=10>]              : Number of iterations");
    addParamsLine("  [--alpha0 <N=0.05>]          : Significance");
    //addParamsLine("  [--alphaF <N=0.005>]         : Final significance");
    //addParamsLine("  [--keepIntermediateVolumes]  : Keep the volume of each iteration");
    addParamsLine("  [--angularSampling <a=5>]    : Angular sampling in degrees for generating the projection gallery");
    addParamsLine("  [--sampling_rate <s=1>]      : Sampling rate in A/px");
    //addParamsLine("  [--maxShift <s=-1>]          : Maximum shift allowed (+-this amount)");
    //addParamsLine("  [--minTilt <t=0>]            : Minimum tilt angle");
    //addParamsLine("  [--maxTilt <t=90>]           : Maximum tilt angle");
    //addParamsLine("  [--useImed]                  : Use Imed for weighting");
    //addParamsLine("  [--strictDirection]          : Images not significant for a direction are also discarded");
    //addParamsLine("  [--angDistance <a=10>]       : Angular distance");

}

void ProgValidationNonTilt::run()
{
    //Clustering Tendency and Cluster Validity Stephen D. Scott

    randomize_random_generator();
    char buffer[400];
    sprintf(buffer, "xmipp_reconstruct_significant -i %s  --initvolumes %s --odir %s --sym  %s --iter 1 --alpha0 %f --angularSampling %f",fnIn.c_str(), fnInit.c_str(),fnDir.c_str(),fnSym.c_str(),alpha0,angularSampling);
    system(buffer);

    MetaData md,mdOut,mdOut2,tempMd2,mdWeight;
    FileName fnMd,fnOut,fnFSC,fnOut2;
    fnMd = fnDir+"/angles_iter01_00.xmd";
    fnOut = fnDir+"/kk.xmd";
    fnOut2 = fnDir+"/kk2.xmd";
    fnFSC = fnDir+"/fsc.xmd";
    size_t nSamplesRandom = 100;

    md.read(fnMd);
    size_t maxNImg;
    size_t sz = md.size();
    md.getValue(MDL_IMAGE_IDX,maxNImg,sz);

    String expression;
    double W,sumW;
    double tempW;
    MDRow row;

    init_progress_bar(maxNImg);
    for (size_t i=0; i<=maxNImg;i++)
    {
    	MetaData tempMd;
        std::vector<double> sum_u(nSamplesRandom);
        std::vector<double> sum_w(nSamplesRandom);
        std::vector<double> H0(nSamplesRandom);
        std::vector<double> H(nSamplesRandom);

        expression = formatString("imageIndex == %lu",i);
        tempMd.importObjects(md, MDExpression(expression));

        if (tempMd.size()==0)
            continue;

        obtainSumU(tempMd,sum_u,H0);
        obtainSumW(tempMd,sum_w,sum_u,H);

        std::sort(H0.begin(),H0.end());
        std::sort(H.begin(),H.end());

        double P = 0.;
        for(int i=0; i<sum_u.size();i++)
        {
            if (H0.at(i) > H.at(i) )
                P += 1.;
        }

        P = (P/nSamplesRandom);
        row.setValue(MDL_IMAGE_IDX,i);
        row.setValue(MDL_WEIGHT,P);
        mdOut.addRow(row);

        sum_u.clear();
		sum_w.clear();
		H0.clear();
		H.clear();
		tempMd.clear();
        progress_bar(i+1);

    }

    mdOut.write(fnOut);

    FileName fnVolume=fnDir+"/volume_projMatch.vol";
    String args=formatString("-i %s -o %s --sym %s --weight -v 0",fnMd.c_str(),fnVolume.c_str(),fnSym.c_str());
    String cmd=(String)"xmipp_reconstruct_fourier "+args;

    if (system(cmd.c_str())==-1)
        REPORT_ERROR(ERR_UNCLASSIFIED,"Cannot open shell");


    // Size of the images
    size_t Xdim, Ydim, Zdim,Ndim;
    getImageSize(fnVolume,Xdim,Ydim,Zdim,Ndim);
    args=formatString("-i %s --mask circular %d -v 0",fnVolume.c_str(),-Xdim/2);
    cmd=(String)"xmipp_transform_mask "+args;
    if (system(cmd.c_str())==-1)
        REPORT_ERROR(ERR_UNCLASSIFIED,"Cannot open shell");

    FileName fnVolumeSig=fnDir+"/volume_iter01_00.vol";
    sprintf(buffer, "xmipp_resolution_fsc --ref %s -i %s -s %f -o %s",fnVolume.c_str(),fnVolumeSig.c_str(),sampling_rate,fnFSC.c_str());
    system(buffer);
}

void ProgValidationNonTilt::obtainSumU(MetaData & tempMd,std::vector<double> & sum_u,std::vector<double> & H0)
{
    //MetaData md,tempMd,mdWeight;
    //FileName fnMd;
    //randomize_random_generator();
    double xRan,yRan,zRan,norm;
    double tilt,rot;
    double sumWRan;
    double xRanArray[tempMd.size()];
    double yRanArray[tempMd.size()];
    double zRanArray[tempMd.size()];
    std::vector<double> weightV;
    double a;

    SymList SL;
    int symmetry, sym_order;
    SL.readSymmetryFile(fnSym.c_str());
    SL.isSymmetryGroup(fnSym.c_str(), symmetry, sym_order);

    double non_reduntant_area_of_sphere = SL.nonRedundantProjectionSphere(symmetry,sym_order);
    double area_of_sphere_no_symmetry = 4.*PI;
    double correction = std::sqrt(non_reduntant_area_of_sphere/area_of_sphere_no_symmetry);

    for (size_t n=0; n<sum_u.size(); n++)
    {
        sumWRan = 0;
        for (size_t nS=0; nS<tempMd.size(); nS++)
        {
             /*
                x = sin(tilt*PI/180)*cos(rot*PI/180);
        		y = sin(tilt*PI/180)*sin(rot*PI/180);
        		z = std::abs(cos(tilt*PI/180));
             */

        	tilt =(double(std::rand())/RAND_MAX)*(PI/2);
        	rot  =(std::rand()-RAND_MAX/2)*(PI/RAND_MAX);
        	xRan = sin(tilt)*cos(rot);
        	yRan = sin(tilt)*sin(rot);
        	zRan = (cos(tilt));

        	//std::cout << tilt << " " << rot << std::endl;
        	//std::cout << xRan << " " << yRan << " " << zRan << " " << std::endl;
            //zRan=(std::rand()-RAND_MAX/2);
            //norm = std::sqrt(xRan*xRan+yRan*yRan+zRan*zRan);
            //xRan = (xRan/norm);
            //yRan = (yRan/norm);
            //zRan = std::abs(zRan/norm);

            xRanArray[nS] = xRan;
            yRanArray[nS] = yRan;
            zRanArray[nS] = zRan;
            tempMd.getColumnValues(MDL_WEIGHT, weightV);

            std::random_shuffle(weightV.begin(), weightV.end());
        }

        sumWRan = 0;
        double WRan, tempWRan, tempW1, tempW2;
        for (size_t nS1=0; nS1<tempMd.size(); nS1++)
        {
            tempWRan = 1e3;
            for (size_t nS2=0; nS2<tempMd.size(); nS2++)
            {
                a = std::abs(std::acos(xRanArray[nS1]*xRanArray[nS2]+yRanArray[nS1]*yRanArray[nS2]+zRanArray[nS1]*zRanArray[nS2]));
                if ( (a<tempWRan) && (a != 0))
                {
                    tempWRan = a;
                    tempW2 = weightV[nS2];
                    tempW1 = weightV[nS1];
                    WRan = a*std::exp(std::abs(tempW1-tempW2))*std::exp(-(tempW1+tempW2));
                    //WRan = a;
                }
            }
            sumWRan += WRan;
        }
        sum_u.at(n)=sumWRan*correction;
    }

    size_t idx = 0;
    while (idx < sum_u.size())
    {
        std::random_shuffle(sum_u.begin(), sum_u.end());

        if(sum_u.at(0) != sum_u.at(1))
        {
            H0[idx] = sum_u.at(0)/(sum_u.at(0)+sum_u.at(1));
            idx += 1;
        }
    }

}

#define _FOR_ALL_OBJECTS_IN_METADATA2(__md) \
        for(MDIterator __iter2(__md); __iter2.hasNext(); __iter2.moveNext())
void ProgValidationNonTilt::obtainSumW(MetaData & tempMd,std::vector<double> & sum_W,std::vector<double> & sum_u,std::vector<double> & H)
{
    double a;
    double rot,tilt,w;
    double x,y,z;
    double xx,yy,zz;
    double w2;
    double tempW;
    double W;
    double sumW;

    sumW = 0;
    FOR_ALL_OBJECTS_IN_METADATA(tempMd)
    {
        tempMd.getValue(MDL_ANGLE_ROT,rot,__iter.objId);
        tempMd.getValue(MDL_ANGLE_TILT,tilt,__iter.objId);
        tempMd.getValue(MDL_WEIGHT,w,__iter.objId);
        x = sin(tilt*PI/180.)*cos(rot*PI/180.);
        y = sin(tilt*PI/180.)*sin(rot*PI/180.);
        z = std::abs(cos(tilt*PI/180.));

        tempW = 1e3;
        _FOR_ALL_OBJECTS_IN_METADATA2(tempMd)
        {
            tempMd.getValue(MDL_ANGLE_ROT,rot,__iter2.objId);
            tempMd.getValue(MDL_ANGLE_TILT,tilt,__iter2.objId);
            tempMd.getValue(MDL_WEIGHT,w2,__iter2.objId);
            xx = sin(tilt*PI/180.)*cos(rot*PI/180.);
            yy = sin(tilt*PI/180.)*sin(rot*PI/180.);
            zz = std::abs(cos(tilt*PI/180.));
            a = std::abs(std::acos(x*xx+y*yy+z*zz));

            if ( (a<tempW) && (a != 0))
            {
                W = a*std::exp(std::abs(w-w2))*std::exp(-(w+w2));
                //W = a;
                tempW = a;
            }
        }
        sumW +=  W;
    }

    size_t idx = 0;
    for (size_t n=0; n<sum_u.size(); n++)
    {
        H[n] = sumW/(sumW+sum_u.at(n));
    }
}

