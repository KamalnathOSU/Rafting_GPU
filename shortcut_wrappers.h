#ifndef SHORTCUT_WRAPPERS_H_INCLUDED
#define SHORTCUT_WRAPPERS_H_INCLUDED

//For looping over the system in real space
#define for_xijk \
for(x_i=0;x_i<L;x_i++) \
	for(x_j=0;x_j<M;x_j++) \
		for(x_k=0;x_k<N;x_k++) \
		{ xijk=x_i*M*N+x_j*N+x_k;
#define efor_xijk }
#define fxijk \
for(x_i=0;x_i<L;x_i++) \
	for(x_j=0;x_j<M;x_j++) \
		for(x_k=0;x_k<N;x_k++) \
		{ 
#define efxijk }
//In the above loop, the order of the looping is important, as we read the input file in the same way


//For looping over the system in k space
#define for_gijk \
for(g_i=0;g_i<L;g_i++) \
	for(g_j=0;g_j<M;g_j++) \
		for(g_k=0;g_k<(N/2+1);g_k++) \
		{ gijk=g_i*M*(N/2+1)+g_j*(N/2+1)+g_k;
#define efor_gijk }
#define fgijk \
for(g_i=0;g_i<L;g_i++) \
	for(g_j=0;g_j<M;g_j++) \
		for(g_k=0;g_k<(N/2+1);g_k++) \
		{ 
#define efgijk }


#endif // SHORTCUT_WRAPPERS_H_INCLUDED
