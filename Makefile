#Updated with the computation of plastic strain
CC=pgcc 

HEADERFILE=constants.h shortcut_wrappers.h global_variables.h Makefile pfr_cufft1.h curand_cuda.h
OBJ=initialization.o pfr_cufft1.o cal_Bpqmatrix.o time_evolution.o io.o cal_kvec.o global_variables.o debug_functions.o initial_conditions.o curand_cuda.o

CFLAGS  = -fast -acc -ta=tesla:cuda9.0 -Minfo=accel -Mcudalib=cufft,curand -Minform=warn -mcmodel medium 

PROJECT_NAME=Ning_raft_full_acc
filename=testrun5.c

%.o: %.c ${HEADERFILE}
	$(CC) -c -o $@ $< $(CFLAGS)

all: ${PROJECT_NAME}

${PROJECT_NAME}: ${filename} $(HEADERFILE) $(OBJ)
	$(CC)  -o $@.exe ${filename} $(OBJ) ${CFLAGS} ${LDFLAGS}

curand_test: curand_test.c curand_cuda.o
	$(CC) $(CFLAGS) -o $@.exe $< curand_cuda.o $(LDFLAGS)

cufft_test: cufft_test.o pfr_cufft1.o
	$(CC) $(CFLAGS) -o $@.exe $< pfr_cufft1.o

test:test.c 
	$(CC) -o $@.exe test.c $(CFLAGS)
clean:
	$(RM) -rf *.o *.exe a.out

cleanall:
	$(RM) -rf *.bin *.dat *.o *.exe a.out *~ output_files/*.bin output_files/*.vtk output_files/*.dat
