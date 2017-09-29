/***************************************************************************
 *
 * Authors:     Carlos Oscar S. Sorzano (coss@cnb.csic.es)
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

#include "volume_texture.h"
#include <data/xmipp_fftw.h>

void ProgVolumeTexture::defineParams()
{
    // Parameters
    addParamsLine(" -i <volume> : Input volume");
    addParamsLine(" -r <volume> : Reference volume");
    addParamsLine(" [--patchSize <n=5>] : Patch size for the correlations");
    mask.defineParams(this,INT_MASK);
}

void ProgVolumeTexture::readParams()
{
	fnIn=getParam("-i");
	fnRef=getParam("-r");
	patchSize=getIntParam("--patchSize");
    mask.allowed_data_types = INT_MASK;
    if (checkParam("--mask"))
        mask.readParams(this);
}

void ProgVolumeTexture::show()
{
    if (verbose)
		std::cout
		<< "Input volume:     " << fnIn   << std::endl
		<< "Reference volume: " << fnRef  << std::endl
		<< "Patch size:       " << patchSize << std::endl
		;
}

//#define DEBUG
void ProgVolumeTexture::run()
{
    show();
    Image<double> V;
    Image<double> R;
    MultidimArray<double> patch;
    int patchSize_2=patchSize/2;
	patch.resize(patchSize,patchSize,patchSize);
    patch.setXmippOrigin();

	V.read(fnIn);
    MultidimArray<double> &mV=V();
    R.read(fnRef);
    MultidimArray<double> &mR=R();

    mask.generate_mask(mV);
	const MultidimArray<int> &mmask = mask.get_binary_mask();

	std::vector< MultidimArray<double> > patchList;
	std::vector< MultidimArray<double> > patchListR;
	int patchListLength=0;
    for (int k=patchSize_2; k<(int)ZSIZE(mV)-patchSize_2; k+=patchSize)
        for (int i=patchSize_2; i<(int)YSIZE(mV)-patchSize_2; i+=patchSize)
            for (int j=patchSize_2; j<(int)XSIZE(mV)-patchSize_2; j+=patchSize)
            	if (mmask(k,i,j)!=0)
		        {
		        	mV.window(patch,
		            		k-patchSize_2,i-patchSize_2,j-patchSize_2,
		            	 	k+patchSize_2-1,i+patchSize_2-1,j+patchSize_2-1);
		            patchList.push_back(patch);

		            mR.window(patch,
		            		k-patchSize_2,i-patchSize_2,j-patchSize_2,
		            	 	k+patchSize_2-1,i+patchSize_2-1,j+patchSize_2-1);
		            patchListR.push_back(patch);

		            patchListLength++;
				}


	std::cout << "Vector length = " << patchListLength  << std::endl
			  << "Patch  xdim  = " << patchList[0].xdim << std::endl
			  << "Patch  ydim  = " << patchList[0].ydim << std::endl
			  << "Patch  zdim  = " << patchList[0].zdim << std::endl
			 ;

	CorrelationAux aux;

	MultidimArray< double> textureCorr, noiseCorr;
	textureCorr.resize(patchSize,patchSize,patchSize);
    textureCorr.setXmippOrigin();
	noiseCorr.resize(patchSize,patchSize,patchSize);
    noiseCorr.setXmippOrigin();

    std::cout << "noiseCorr size = " <<  noiseCorr.xdim << "x" 
     		  << noiseCorr.ydim << "x" << noiseCorr.zdim << std::endl
			  << "textureCorr size = " <<  textureCorr.xdim << "x" 
     		  << textureCorr.ydim << "x" << textureCorr.zdim << std::endl  ;

    int n;
	double cumTextureCorr=0, cumNoiseCorr=0;
	for (int i=0; i<patchListLength ; i++)
		for (int j=0; j<patchListLength ; j++)
		{
			n = i+j;
			if(n>=patchListLength) n-=patchListLength;
			correlation_matrix(patchList[i],patchList[n],noiseCorr,aux,false);
			correlation_matrix(patchList[i],patchListR[n],textureCorr,aux,false);
			cumTextureCorr += textureCorr.sum();
			cumNoiseCorr += noiseCorr.sum();
		}
	
	std::cout << "Cumulative TextureCorrelation        -> " << cumTextureCorr << std::endl
			  << "Cumulative NoiseCorrelation          -> " << cumNoiseCorr << std::endl 
			  << "Cumulative Texture/Noise Correlation -> " << cumTextureCorr/cumNoiseCorr
			  << std::endl ;

}
#undef DEBUG
