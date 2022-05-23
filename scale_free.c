/**
 * @file main.c
 * @author Rilder Pires
 * @version 1.0
 *
 * @section LICENSE
 *
 * Barabasi Albert Algorithm
 *
 * Copyright (C) 2013  Rilder Pires
 * All rights reserved.
 *
 * Copyright (C) 1997 - 2002, Makoto Matsumoto and Takuji Nishimura,
 * All rights reserved.
 *
 * Copyright (C) 2005, Mutsuo Saito,
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 *  1. Redistributions of source code must retain the above copyright
 *     notice, this list of conditions and the following disclaimer.
 *
 *  2. Redistributions in binary form must reproduce the above copyright
 *     notice, this list of conditions and the following disclaimer in the
 *     documentation and/or other materials provided with the distribution.
 *
 *  3. The names of its contributors may not be used to endorse or promote
 *     products derived from this software without specific prior written
 *     permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * @section DESCRIPTION
 *
 * This is an C implementation of Barabasi Albert algorithm coded by Rilder
 * Pires in 2013/08/30. The original idea of this algorithm was proposed by
 * A. Barabási and R. Albert in 1999. It was decribed into a paper called
 * "Emergence of Scaling in Random Networks".
 *
 * Any feedback is very welcome.
 * e-mail: rilder@fisica.ufc.br
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "mtwister.h"

int main(){
	int i,j,ii,cv,ok;
	double r;
	int m0=5;
	int m=5;
	int t=150000;

	int edges[m*t];
	int k[m0+t];
	//int *edges;
	//int *k;
	//int *dist;
	int kmax=0;

	/* Set the seed of Merssene Twister */
	MTGenSeed(764763469);

	/* Alloc memory to keep track of the vertices that are linked with each vertex */
	//edges=(int*)calloc(m*t,sizeof(int));

	/* Alloc memory to keep track of the conectivity of each vertex */
	//k=(int*)calloc(m0+t,sizeof(int));

	/* Edges of the network */
	/* edges[i*m+j] (0 <= j < m) are the j vertices that are linked to i. */
	for(j=0;j<m*t;j++) edges[j]=-1;

	/* This statement deal with the inclusion of the first vertex of the network */
	for(j=0;j<m;j++){
		edges[j]=j;
		k[j]++;
		k[m0+0]++;
		printf("%d %d %d %d\n",j, edges[j],k[j],k[m0+0]);
	};
	kmax=k[m0+0];

	/* This statement deal with the inclusion of the other vertices of the network */
	for(i=1;i<t;i++){
		/* Put m edges to the i vertex. */
		for(ii=0;ii<m;ii++){
			ok=0;
			/* Sort a vertex until it satisfy our conditions (Regection Method).*/
			while(ok==0){
				cv=MTrand_int31()%(m0+i);
				r=MTrandReal3()*((double)kmax);
				/* See if the choosen point is below the probability curve. */
				if(r<((double)k[cv])){
					ok=1;
					for(j=0;j<ii;j++){
						if(edges[i*m+j]==cv) ok=0;
					};
				};
			};
			/* The vertex was chosen, so now is time to include it to our list. */
			edges[i*m+ii]=cv;
			k[cv]++;
			k[m0+i]++;
			/* Here we store the value of the maximum of the probability distribution. */
			if(k[cv]>kmax) kmax=k[cv];
			if(k[m0+i]>kmax) kmax=k[m0+i];
		};
	};

	/* Deal with statistics stuf (in the future it's a good choice make this in
	other program) */
	//dist=(int*)calloc(kmax+1,sizeof(int));
	int dist[kmax+1];

	ii=0;
	for(i=0;i<m0+t;i++){
		dist[k[i]]++;
		ii+=k[i];
		if(k[i]==-1) printf("%d\n",k[i]);
	};

	j=0;
	for(i=0;i<kmax+1;i++){
		printf("%d %d\n",i,dist[i]);
		j+=dist[i];
	};

	//free(edges);
	//free(dist);
	//free(k);
	return 0;
}
