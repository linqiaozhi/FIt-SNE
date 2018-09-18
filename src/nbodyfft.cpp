#include "winlibs/stdafx.h"
#include "parallel_for.h"

#include "nbodyfft.h"


using namespace std;

long int diff(timespec start, timespec end)
{
    long int tv_nsec;
    tv_nsec = 1000000000*(end.tv_sec-start.tv_sec)+end.tv_nsec-start.tv_nsec;
    return tv_nsec;
}

int precompute2( double xmax, double xmin, double ymax, double ymin,  int nlat,
		int nterms, kerneltype2d ker, double * band,double * boxl,
		double * boxr,  double * prods, double * xpts, double *
		xptsall, double *yptsall, int * irearr, fftw_complex * zkvalf )
{

	/*
	 * Set up the boxes
	 */
	int nboxes = nlat*nlat;
	double boxwidth = (xmax - xmin)/(double) nlat;

	//Left and right bounds of each box
	int nn = nboxes*nterms*nterms;
	if (boxl == NULL) {
		printf("Malloc failed\n");
		exit(-1);
	}

	int ii=0;
	//printf("boxwidth: %lf, nlat %d, xmax %f xmin %f\n", boxwidth, nlat, xmax, xmin);
	for (int i = 0; i < nlat; i++){
		for (int j = 0; j < nlat; j++){
			boxl[0*nboxes + ii] = j*boxwidth + xmin;
			boxr[0*nboxes + ii] = (j+1)*(boxwidth) + xmin;

			boxl[1*nboxes + ii] = i*boxwidth + ymin;
			boxr[1*nboxes + ii] = (i+1)*(boxwidth) + ymin;
			//printf("box %d, %lf to %lf and %lf %lf \n", i, boxl[i], boxr[i] , boxl[nboxes+ii], boxr[nboxes+ii]);
			ii++;
		}
	}

	//Coordinates of each (equispaced) interpolation node for a single box
	double h = 2/(double)nterms;

	xpts[0] = -1 + h/2.0;
	for (int i=1; i< nterms; i++){
		xpts[i] = xpts[0] + (i)*h;
	}
	/*
	 * Interpolate kernel using lagrange polynomials
	 */

	//Get the denominators for the lagrange polynomials (this is getprods())
	for (int i= 0; i< nterms; i ++ ){
		prods[i] = 1;
		for (int j=0;j < nterms; j ++ ){
			if (i != j) {
				prods[i] = prods[i] * (xpts[i] - xpts[j]);
			}
		}
                prods[i] = 1/prods[i];
//		printf("Prods[%d] xpts[%d] = %lf, %lf\n", i,i, prods[i], xpts[i]);
	}


	//Coordinates of each (equispaced) interpolation node for all boxes
	int nfourh = nterms*nlat;
	int nfour = 2*nterms*nlat;
	h = h*boxwidth/2;
	double xstart = xmin + h/2;
	double ystart = ymin + h/2;


	ii = 0;
	for (int i= 0; i< nfourh; i ++ ){
		for (int j= 0; j< nfourh; j ++ ){
			xptsall[ii] = xstart + (i)*h;
			yptsall[ii] = ystart + (j)*h;
			//printf("%d xptsall %lf, %lf \n", ii, xptsall[ii],  yptsall[ii]);
			ii++;
		}
	}

	ii = 0;
	for (int ilat = 0; ilat < nlat; ilat++){
		for (int jlat = 0; jlat < nlat; jlat++){
			for (int i = 0; i < nterms; i++){
				for (int j = 0; j < nterms; j++){
					int iy = (ilat)*nterms + j;
					int ix = (jlat)*nterms + i;
					int ipt = (ix)*nlat*nterms + iy;
					irearr[ii] = ipt;
					//printf("irearr[%d]=%d\n", ii, ipt);
					ii++;
				}
			}
		}
	}


	//Kernel evaluated at interpolation nodes. Make it circulant
	double * zkvals = (double *) calloc(nfour*nfour,sizeof(double));
	ii = 0;
	for (int i = 0; i< nfourh; i++){
		for (int j = 0; j< nfourh; j++){
			double tmp = ker(xptsall[0],yptsall[0], xptsall[ii],yptsall[ii], band[i], band[i]);
			//printf("tmp[%d,%d] %f, %f, %f,%f = %f\n", i,j,xptsall[0],yptsall[0], xptsall[ii],yptsall[ii],tmp);
			zkvals[(i+nfourh)*nfour + (j+nfourh)] = tmp;
			zkvals[(nfourh-i)*nfour + (j+nfourh)] = tmp;
			zkvals[(i+nfourh)*nfour + (nfourh -j )] = tmp;
			zkvals[(nfourh-i)*nfour + (nfourh-j)] = tmp;

			ii++;
		}
	}
	for (int i = 20; i< 30; i++){
		for (int j = 30; j< 31; j++){

			//printf("zkvals[%d,%d] = %f\n", i,j,zkvals[i*nfour + j]);
		}
	}

	//FFT of the kernel
	double * zkvali = (double*) fftw_malloc(sizeof(double) * nfour*nfour);
	fftw_plan p;
	p = fftw_plan_dft_r2c_2d(nfour,nfour, zkvali, zkvalf, FFTW_ESTIMATE);
	for (int i =0; i< nfour*nfour; i++){
		zkvali[i] =  0;
	}

	for (int i = 0; i< nfour*nfour; i++){
		zkvali[i] =  zkvals[i];
	}
	fftw_execute(p);

	/*
	for (int i = 20; i< 30; i++){
		for (int j = 30; j< 31; j++){
			printf("zkvalsf[%d,%d] = %f\n", i,j,zkvalf3[i*nfour + j][0]/(nfour*nfour));
		}
	}
	*/

	fftw_destroy_plan(p);
	fftw_free(zkvali);

	free(zkvals);
	//The rest of this should be in a separate function...

	return 1;
}


int nbodyfft2(int n, int ndim, double* xs, double *ys, double * charges, int
		nlat, int nterms,double *boxl, double *boxr,  double * prods,
		double * xpts, double * xptsall, double *yptsall,int* irearr,
		fftw_complex * zkvalf, double * outpot,unsigned int nthreads ){

    if (nthreads == 0) {
        nthreads = std::thread::hardware_concurrency();
    }

        struct timespec start10, end10, start20, end20,start30, end30;
        clock_gettime(CLOCK_MONOTONIC, &start10);
	int nboxes = nlat*nlat;
	int nn = nboxes*nterms*nterms;

	double rmin = boxl[0];
	double boxwidth = boxr[0] - boxl[0];

	int * boxcount =(int*) calloc((nboxes +1), sizeof(int));
	int * boxcounti =(int*) calloc(nboxes +1, sizeof(int));
	int * boxsort =(int*) malloc(n* sizeof(int));
	if (boxcount == NULL) {
		printf("Malloc failed\n");
		exit(-1);
	}


	//Initialize the charges and the locations
	for (int i=0; i< n; i++){
		int ixs = (xs[i] - rmin)/boxwidth;
		int iys = (ys[i] - rmin)/boxwidth;
		//printf("%lf,%lf\n", rmin, boxwidth);
		//printf("%d,%d,%d\n", ixs, iys, nboxes);
		if (ixs >= nlat) {
			ixs = nlat - 1;
		}
		if (ixs < 0) {
			ixs = 0;
		}
		if (iys >= nlat) {
			iys = nlat - 1;
		}
		if (iys < 0) {
			iys = 0;
		}

		int icol =  ixs;
		int irow =  iys;

		int iadr = (irow)*nlat + icol;

		boxcount[iadr] += 1;
		boxsort[i] = iadr;
		//printf("%d: %f,%f, in box %d, x: %.2f - %.2f,y: %.2f - %.2f which has %d\n", i, xs[i],ys[i],  iadr, boxl[0*nboxes + iadr], boxr[0*nboxes + iadr],boxl[1*nboxes + iadr], boxr[1*nboxes + iadr], boxcount[iadr]);
	}


	int * iarr =(int*) malloc(n* sizeof(int));
	double * chargessort =(double*) malloc(ndim*n* sizeof(double));

	//boxsort[i] = ibox: the box for the ith point
	//Set the offset of each box
	int * boxoffset =(int*) malloc((nboxes+1)* sizeof(int));
	boxoffset[0] = 0;
	for (int ibox = 1; ibox<nboxes +1; ibox++){
		boxoffset[ibox] = boxoffset[ibox-1] + boxcount[ibox-1];
	}

	for (int ibox = 0; ibox<nboxes; ibox++){
		boxcounti[ibox] = 0;
	}

	//The number of points in each box (so far)
	for (int i=0; i<n; i++){
		int iadr =boxsort[i];
		int indx = boxoffset[iadr] + boxcounti[iadr];
		iarr[indx] = i;
		boxcounti[iadr] = boxcounti[iadr] +1;
		//printf("%f, %f , iarr[%d] = %d, iadr = %d, boxoffset[iadr] = %d, boxcounti[iadr] = %d\n", xs[iarr[indx]],ys[iarr[indx]],  indx,i,iadr, boxoffset[iadr], boxcounti[iadr]);
	}

	   for (int i=0; i<n; i++){
	   int ibox = boxsort[iarr[i]];
	   //printf("%d (%d): %f, in box %d, %.2f - %.2f which has %d\n", i,iarr[i], locs[iarr[i]], ibox, boxl[ibox], boxr[ibox],boxcount[ibox]);
	   }

	//Locsort

	//FILE *f = fopen("iarr.txt", "w");
        //Parallelize. separate charges
	double * xsort =(double*) malloc(n* sizeof(double));
	double * ysort =(double*) malloc(n* sizeof(double));



    	PARALLEL_FOR(nthreads,n, {
		xsort[loop_i] = xs[iarr[loop_i]];
		ysort[loop_i] = ys[iarr[loop_i]];
		});
		//fprintf(f, "%d,", iarr[i]);
		for (int idim=0; idim<ndim; idim++){
      PARALLEL_FOR(nthreads,n,{
			chargessort[idim*n+loop_i] = charges[idim*n +iarr[loop_i]];
	    });
		}

	for (int i=0; i<10; i++){
//		printf("Charge %d at %f,%f, sorted %f,%f: %f, sorted; %f\n", i, xs[i], ys[i], xsort[i],xsort[i], charges[i], chargessort[0*n+i]);
	}



	//tlocssort is the translated locations
        //396 in fortran code. not worth parallelizing. Pull out the boxl
	double * xsp = (double *) malloc(n*sizeof(double));
	double * ysp = (double *) malloc(n*sizeof(double));
	double bsizeinv = 2/(boxr[1]-boxl[1]);
	for (int ibox=0; ibox<nboxes;ibox++){
			double xmin = boxl[ibox];
			double xmax = boxr[ibox];
			double ymin = boxl[nboxes+ ibox];
			double ymax = boxr[nboxes+ ibox];
		for (int i=boxoffset[ibox]; i<boxoffset[ibox+1];i++){
			xsp[i] = (xsort[i] - xmin)*bsizeinv - 1;
			ysp[i] = (ysort[i] - ymin)*bsizeinv - 1;
			//printf("i %d  %f,%f  tlocssort[i] %f,%f  xmin %f xmax %f, ymin %f ymax %f\n", i, xsort[i],ysort[i],  xsp[i],ysp[i], xmin, xmax, ymin, ymax);
		}
	}

        clock_gettime(CLOCK_MONOTONIC, &end10);
        printf("Initializing and sorting (%d threads): %.2lf ms\n",nthreads, (diff(start10,end10))/(double)1E6);

	//Get the L_j vals
        clock_gettime(CLOCK_MONOTONIC, &start20);

        //parallelize
	double * ydiff = (double*) malloc(n*nterms*sizeof(double));
	double * yprods = (double*) malloc(n*nterms*sizeof(double));

	    PARALLEL_FOR(nthreads,n,{
		yprods[loop_i] = 1;
		for (int i =0; i < nterms; i++) {
			ydiff[loop_i + i*n] = xsp[loop_i] - xpts[i];
			yprods[loop_i] = yprods[loop_i]*ydiff[loop_i+i*n];
		}
	    });

        //compute inverse of prods in line 66 of my code and just multiply
        //Parallelize
	double * svalsx = (double*) malloc(n*nterms*sizeof(double));


	PARALLEL_FOR(nthreads,n, {
	for (int i =0; i < nterms; i++) {
		if ( fabs(ydiff[loop_i+i*n]) >= 1e-6) {
			svalsx[loop_i+i*n] = yprods[loop_i]*prods[i]/ydiff[loop_i+i*n];
		}
		if ( fabs(ydiff[loop_i+i*n]) < 1e-6) {
			svalsx[loop_i+i*n] =  prods[i];
			for (int k =0; k < nterms; k++) {
				if(i != k) {
					svalsx[loop_i+i*n] = svalsx[loop_i+i*n] *ydiff[loop_i+k*n];
				}
			}
		}
	}
    });

	//L_j for y
	for (int j =0; j < n; j++) {
		yprods[j] = 1;
		for (int i =0; i < nterms; i++) {
			ydiff[j + i*n] = ysp[j] - xpts[i];
			yprods[j] = yprods[j]*ydiff[j+i*n];
		}
	}


        //Make prods be 1/prods
	double * svalsy = (double*) malloc(n*nterms*sizeof(double));
                    PARALLEL_FOR(nthreads,n, {
                                    for (int i =0; i < nterms; i++) {
                                            if ( fabs(ydiff[loop_i+i*n]) >= 1e-6) {
                                                    svalsy[loop_i+i*n] = yprods[loop_i]*prods[i]/ydiff[loop_i+i*n];
                                                            //printf("svals[%d] = %lf\n", loop_i+i*n, svals[loop_i+i*n]);
                                            }
                                            if ( fabs(ydiff[loop_i+i*n]) < 1e-6) {
                                                    svalsy[loop_i+i*n] =  prods[i];
                                                    for (int k =0; k < nterms; k++) {
                                                            if(i != k) {
                                                                    svalsy[loop_i+i*n] = svalsy[loop_i+i*n] *ydiff[loop_i+k*n];
                                                            }
                                                    }
                                            }
                                    }
                    })
    	    clock_gettime(CLOCK_MONOTONIC, &end20);
        printf("Getting svals (%d threads): %.2lf ms\n", nthreads,(diff(start20,end20))/(double)1E6);

	//Compute mpol, which is interpolation


        struct timespec start, end, start2, end2,start3, end3;
        clock_gettime(CLOCK_MONOTONIC, &start);

	double * mpol = (double *) calloc(sizeof(double), nn*ndim);
	double * loc = (double *) calloc(sizeof(double), nn*ndim);



	PARALLEL_FOR(nthreads,nboxes,
                            {
                                    int istart = loop_i*nterms*nterms;
                                    for (int impx=0; impx<nterms; impx++){
                                            for (int impy=0; impy<nterms; impy++){
                                                int istartii = istart+(impx)*nterms + impy;
                                                for (int idim=0;idim<ndim; idim++){
                                                        double temp2 = 0;
                                                        for (int i=boxoffset[loop_i]; i<boxoffset[loop_i+1];i++){
                                                                double temp = svalsy[impy*n + i]*svalsx[impx*n + i];
                                                                temp2 += temp*chargessort[idim*n+i];
                                                            }
                                                            mpol[idim + (istartii)*ndim] += temp2;
                                                    }
                                            }
                                    }
                            });

    clock_gettime(CLOCK_MONOTONIC, &end);
    printf("Step 1 loop (%d threads): %.2lf ms\n", nthreads, (diff(start,end))/(double)1E6);


	//Mpol to loc

	int nfourh = nterms*nlat;
	int nfour = 2*nterms*nlat;

    clock_gettime(CLOCK_MONOTONIC, &start2);
		fftw_plan p,p2;
		//FFT of zmpol
		//printf("doing %d, by %d\n", nfour, nfour);
		//fftw_complex * zmpoli = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nfour*nfour);
		//fftw_complex * zmpolf = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nfour*nfour);
		//p = fftw_plan_dft_2d(nfour,nfour, zmpoli, zmpolf, FFTW_FORWARD, FFTW_ESTIMATE);
		//p2 = fftw_plan_dft_2d(nfour,nfour, zmpolf, zmpolf, FFTW_BACKWARD, FFTW_ESTIMATE);

		//fftw_import_wisdom_from_string(wisdom_string);

		double * zmpoli = (double*) fftw_malloc(sizeof(fftw_complex) * nfour*nfour);
		fftw_complex * zmpolf = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nfour*(nfour/2 +1));
		double * zmpolfo = (double*) fftw_malloc(sizeof(fftw_complex) * nfour*nfour);
		p = fftw_plan_dft_r2c_2d(nfour,nfour, zmpoli, zmpolf, FFTW_ESTIMATE );
		p2 = fftw_plan_dft_c2r_2d(nfour,nfour, zmpolf, zmpolfo, FFTW_ESTIMATE);
		//wisdom_string = fftw_export_wisdom_to_string();
		double * mpolsort = (double *) calloc(sizeof(double), nn);
		double * zmpol = (double *) calloc(sizeof(double), nfour*nfour);

	for (int idim =0; idim< ndim; idim++){
		for (int i =0; i< nn; i++){
			mpolsort[irearr[i]] = mpol[i*ndim+idim];
		}
		for (int j =0; j < 100; j++) {
			//	printf("mpolsort[%d]=%lf\n",j, mpolsort[j]);
		}

		for (int i =0; i< nfourh; i++){
			for (int j =0; j< nfourh; j++){
				int ii =  i*nfourh +j;
				zmpol[i*nfour + j] = mpolsort[ii];
				//printf("zmpol[%d,%d] = %lf\n", i,j , mpolsort[ii]);
			}
		}



		for (int i =0; i< nfour*nfour; i++){
			zmpoli[i] =  0;
		}
		for (int i = 0; i< nfour*nfour; i++){
			zmpoli[i] =  zmpol[i];
		}

		fftw_execute(p);

		//Take hadamard product

		for (int i =0; i< nfour*(nfour/2 +1); i++){
			//(x_ + y_*i) * (u_ + v_*i) = (x*u - y*v) + (x*v+y*u)i
			double x_ = zmpolf[i][0];
			double y_ = zmpolf[i][1];
			double u_ = zkvalf[i][0];
			double v_ = zkvalf[i][1];
			zmpolf[i][0] = (x_*u_ - y_*v_);
			zmpolf[i][1] = (x_*v_ + y_*u_);
					//printf("(%lf + %lfi) (%lf + %lfi) = %lf + %lfi\n", x_, y_, u_,v_, zmpolf[i][0], zmpolf[i][1]);
		}
		for (int i=0;i<5; i++){
			for (int j=0;j<5; j++){
				//printf("before zmpolf[%d,%d : %lf\n ", i,j, zmpolf[i*nfour+j][0]);
			}
		}

		//Inverse it!
		fftw_execute(p2);
		for (int i=0;i<nfour; i++){
			for (int j=0;j<nfour; j++){
				//printf("zmpolf[%d,%d : %lf\n ", i,j, zmpolf[i*nfour+j][0]/(double)(nfour*nfour));
				//printf("%0.4e,  ", zmpolf[i*nfour+j][0]/(double)(nfour*nfour));
			}
			//printf("\n%d", i+1);
		//	printf("\n", i+1);
		}

		for (int i=0;i<nfourh; i++){
			for (int j=0;j<nfourh; j++){
				int ii = i*nfourh+j;
				int rowval = (nfourh+i);
				int colval = (nfourh+j);

				mpolsort[ii] = zmpolfo[rowval*nfour+colval]/(double)(nfour*nfour);
				//printf("%d row: %d col %d is %lf\n", ii, rowval, colval, mpolsort[ii]);
			}
		}
		for (int i =0; i< nfourh*nfourh; i++){
			loc[i*ndim + idim] = mpolsort[irearr[i]];
			//printf("loc[%d]: %lf\n", i, loc[i]);
		}

	}
		fftw_free(zmpoli);
		fftw_free(zmpolf);
		fftw_free(zmpolfo);
		fftw_destroy_plan(p);
		fftw_destroy_plan(p2);
		free(zmpol);
		free(mpolsort);


    clock_gettime(CLOCK_MONOTONIC, &end2);
    printf("Step 2 FFT (not threaded): %.2lf ms\n", (diff(start2,end2))/(double)1E6);

    //boxes, j, then l, i, then ndim, 3 or 5 ms.
    clock_gettime(CLOCK_MONOTONIC, &start3);
	double * pot = (double *) calloc(n*ndim,sizeof(double));
			PARALLEL_FOR(nthreads,nboxes,
                            {
                                int istart = loop_i*nterms*nterms;
                                for (int i=boxoffset[loop_i]; i<boxoffset[loop_i+1];i++){
                                        
                                        for (int idim=0;idim<ndim; idim++){
                                            double temp2 =0;
                                            for (int j=0; j<nterms; j++){
                                                for (int l=0; l<nterms; l++){
                                                        int tempii = (istart+j*nterms + l)*ndim;
                                                        double temp = svalsx[j*n + i]*svalsy[l*n+i];
                                                        temp2  += temp*loc[tempii+idim];
                                                    }
                                            }
                                            pot[i*ndim +idim] += temp2;
                                    }
                                }
                            });





	for (int i=0; i<n;i++){
		//printf("pot[%d]= %lf\n", i, pot[i]);
	}

//parallelize
	PARALLEL_FOR(nthreads,n,{
		for (int j=0; j<ndim;j++){
			outpot[iarr[loop_i]*ndim+j] = pot[loop_i*ndim+j];
		}
	});
    clock_gettime(CLOCK_MONOTONIC, &end3);
    printf("Step 3 (with unsort) (%d threads): %.2lf ms\n", nthreads, (diff(start3,end3))/(double)1E6);
	free(boxcount); free(boxcounti); free(boxsort); free(iarr); free(chargessort);
	free(boxoffset); free(ydiff); free(yprods);free(svalsx);free(svalsy);free(mpol); free(loc);
	free(pot);
	free(xsort); free(ysort); free(xsp); free(ysp);

	return 1;

}


void precompute(double y_min, double y_max, int n_boxes, int n_interpolation_points, kernel_type kernel,
                double *box_lower_bounds, double *box_upper_bounds, double *y_tilde_spacing, double *y_tilde,
                complex<double> *fft_kernel_vector) {
    /*
     * Set up the boxes
     */
    double box_width = (y_max - y_min) / (double) n_boxes;
    // Compute the left and right bounds of each box
    for (int box_idx = 0; box_idx < n_boxes; box_idx++) {
        box_lower_bounds[box_idx] = box_idx * box_width + y_min;
        box_upper_bounds[box_idx] = (box_idx + 1) * box_width + y_min;
    }

    int total_interpolation_points = n_interpolation_points * n_boxes;
    // Coordinates of each equispaced interpolation point for a single box. This equally spaces them between [0, 1]
    // with equal space between the points and half that space between the boundary point and the closest boundary point
    // e.g. [0.1, 0.3, 0.5, 0.7, 0.9] with spacings [0.1, 0.2, 0.2, 0.2, 0.2, 0.1], respectively. This ensures that the
    // nodes will still be equispaced across box boundaries
    double h = 1 / (double) n_interpolation_points;
    y_tilde_spacing[0] = h / 2;
    for (int i = 1; i < n_interpolation_points; i++) {
        y_tilde_spacing[i] = y_tilde_spacing[i - 1] + h;
    }

    // Coordinates of all the equispaced interpolation points
    h = h * box_width;
    y_tilde[0] = y_min + h / 2;
    for (int i = 1; i < total_interpolation_points; i++) {
        y_tilde[i] = y_tilde[i - 1] + h;
    }

    /*
     * Evaluate the kernel at the interpolation nodes and form the embedded generating kernel vector for a circulant
     * matrix
     */
    auto *kernel_vector = new complex<double>[2 * total_interpolation_points]();
    // Compute the generating vector x between points K(y_i, y_j) where i = 0, j = 0:N-1
    // [0 0 0 0 0 5 4 3 2 1] for linear kernel
    // This evaluates the Cauchy kernel centered on y_tilde[0] to all the other points
    for (int i = 0; i < total_interpolation_points; i++) {
        kernel_vector[total_interpolation_points + i].real(kernel(y_tilde[0], y_tilde[i]));
    }
    // This part symmetrizes the vector, this embeds the Toeplitz generating vector into the circulant generating vector
    // but also has the nice property of symmetrizing the Cauchy kernel, which is probably planned
    // [0 1 2 3 4 5 4 3 2 1] for linear kernel
    for (int i = 1; i < total_interpolation_points; i++) {
        kernel_vector[i].real(kernel_vector[2 * total_interpolation_points - i].real());
    }

    // Precompute the FFT of the kernel generating vector
    fftw_plan p = fftw_plan_dft_1d(2 * total_interpolation_points, reinterpret_cast<fftw_complex *>(kernel_vector),
                                   reinterpret_cast<fftw_complex *>(fft_kernel_vector), FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(p);
    fftw_destroy_plan(p);

    delete[] kernel_vector;
}


void interpolate(int n_interpolation_points, int N, const double *y_in_box, const double *y_tilde_spacings,
                 double *interpolated_values) {
    // The denominators are the same across the interpolants, so we only need to compute them once
    auto *denominator = new double[n_interpolation_points];
    for (int i = 0; i < n_interpolation_points; i++) {
        denominator[i] = 1;
        for (int j = 0; j < n_interpolation_points; j++) {
            if (i != j) {
                denominator[i] *= y_tilde_spacings[i] - y_tilde_spacings[j];
            }
        }
    }
    // Compute the numerators and the interpolant value
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < n_interpolation_points; j++) {
            interpolated_values[j * N + i] = 1;
            for (int k = 0; k < n_interpolation_points; k++) {
                if (j != k) {
                    interpolated_values[j * N + i] *= y_in_box[i] - y_tilde_spacings[k];
                }
            }
            interpolated_values[j * N + i] /= denominator[j];
        }
    }

    delete[] denominator;
}


void nbodyfft(int N, int n_terms, double *Y, double *chargesQij, int n_boxes, int n_interpolation_points,
              double *box_lower_bounds, double *box_upper_bounds, double *y_tilde_spacings, double *y_tilde,
              complex<double> *fft_kernel_vector, double *potentialsQij) {
    int total_interpolation_points = n_interpolation_points * n_boxes;

    double coord_min = box_lower_bounds[0];
    double box_width = box_upper_bounds[0] - box_lower_bounds[0];

    // Determine which box each point belongs to
    auto *point_box_idx = new int[N];
    for (int i = 0; i < N; i++) {
        auto box_idx = static_cast<int>((Y[i] - coord_min) / box_width);
        // The right most point maps directly into `n_boxes`, while it should belong to the last box
        if (box_idx >= n_boxes) {
            box_idx = n_boxes - 1;
        }
        point_box_idx[i] = box_idx;
    }

    // Compute the relative position of each point in its box in the interval [0, 1]
    auto *y_in_box = new double[N];
    for (int i = 0; i < N; i++) {
        int box_idx = point_box_idx[i];
        double box_min = box_lower_bounds[box_idx];
        y_in_box[i] = (Y[i] - box_min) / box_width;
    }

    /*
     * Step 1: Interpolate kernel using Lagrange polynomials and compute the w coefficients
     */
    // Compute the interpolated values at each real point with each Lagrange polynomial
    auto *interpolated_values = new double[n_interpolation_points * N];
    interpolate(n_interpolation_points, N, y_in_box, y_tilde_spacings, interpolated_values);

    auto *w_coefficients = new double[total_interpolation_points * n_terms]();
    for (int i = 0; i < N; i++) {
        int box_idx = point_box_idx[i] * n_interpolation_points;
        for (int interp_idx = 0; interp_idx < n_interpolation_points; interp_idx++) {
            for (int d = 0; d < n_terms; d++) {
                w_coefficients[(box_idx + interp_idx) * n_terms + d] +=
                        interpolated_values[interp_idx * N + i] * chargesQij[i * n_terms + d];
            }
        }
    }

    // `embedded_w_coefficients` is just a vector of zeros prepended to `w_coefficients`, this (probably) matches the
    // dimensions of the kernel matrix K and since we embedded the generating vector by prepending values, we have to do
    // the same here
    auto *embedded_w_coefficients = new double[2 * total_interpolation_points * n_terms]();
    for (int i = 0; i < total_interpolation_points; i++) {
        for (int d = 0; d < n_terms; d++) {
            embedded_w_coefficients[(total_interpolation_points + i) * n_terms + d] = w_coefficients[i * n_terms + d];
        }
    }

    /*
     * Step 2: Compute the values v_{m, n} at the equispaced nodes, multiply the kernel matrix with the coefficients w
     */
    auto *fft_w_coefficients = new complex<double>[2 * total_interpolation_points];
    auto *y_tilde_values = new double[total_interpolation_points * n_terms]();

    fftw_plan plan_dft, plan_idft;
    plan_dft = fftw_plan_dft_1d(2 * total_interpolation_points, reinterpret_cast<fftw_complex *>(fft_w_coefficients),
                                reinterpret_cast<fftw_complex *>(fft_w_coefficients), FFTW_FORWARD, FFTW_ESTIMATE);
    plan_idft = fftw_plan_dft_1d(2 * total_interpolation_points, reinterpret_cast<fftw_complex *>(fft_w_coefficients),
                                 reinterpret_cast<fftw_complex *>(fft_w_coefficients), FFTW_BACKWARD, FFTW_ESTIMATE);

    for (int d = 0; d < n_terms; d++) {
        for (int i = 0; i < 2 * total_interpolation_points; i++) {
            fft_w_coefficients[i].real(embedded_w_coefficients[i * n_terms + d]);
        }
        fftw_execute(plan_dft);

        // Take the Hadamard product of two complex vectors
        for (int i = 0; i < 2 * total_interpolation_points; i++) {
            double x_ = fft_w_coefficients[i].real();
            double y_ = fft_w_coefficients[i].imag();
            double u_ = fft_kernel_vector[i].real();
            double v_ = fft_kernel_vector[i].imag();
            fft_w_coefficients[i].real(x_ * u_ - y_ * v_);
            fft_w_coefficients[i].imag(x_ * v_ + y_ * u_);
        }

        // Invert the computed values at the interpolated nodes, unfortunate naming but it's better to do IDFT inplace
        fftw_execute(plan_idft);

        for (int i = 0; i < total_interpolation_points; i++) {
            // FFTW doesn't perform IDFT normalization, so we have to do it ourselves. This is done by multiplying the
            // result with the number of points in the input
            y_tilde_values[i * n_terms + d] = fft_w_coefficients[i].real() / (total_interpolation_points * 2.0);
        }
    }

    fftw_destroy_plan(plan_dft);
    fftw_destroy_plan(plan_idft);
    delete[] fft_w_coefficients;

    /*
     * Step 3: Compute the potentials \tilde{\phi}
     */
    for (int i = 0; i < N; i++) {
        int box_idx = point_box_idx[i] * n_interpolation_points;
        for (int j = 0; j < n_interpolation_points; j++) {
            for (int d = 0; d < n_terms; d++) {
                potentialsQij[i * n_terms + d] +=
                        interpolated_values[j * N + i] * y_tilde_values[(box_idx + j) * n_terms + d];
            }
        }
    }

    delete[] point_box_idx;
    delete[] y_in_box;
    delete[] interpolated_values;
    delete[] w_coefficients;
    delete[] y_tilde_values;
    delete[] embedded_w_coefficients;
}
