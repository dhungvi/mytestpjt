CC = g++
DEBUG = -g
OPTS = ${DEBUG} -Wall 

O_FILES = Alignment.o AlignmentAlgorithm.o CalculatedOvlDistance.o D2.o ESTAssembly.o \
	  Graph.o InclusionNodes.o OvlDistance.o Prim.o Reconstruction.o SingleBase.o SixTuple.o \
	  SixTuplesGeneration.o 


Main: Main.cpp ${O_FILES}
	${CC} -o Main Main.cpp ${O_FILES} ${OPTS}

Alignment.o: Alignment.cpp Alignment.h Param.h AlignmentAlgorithm.h
	${CC} -c Alignment.cpp ${OPTS}

AlignmentAlgorithm.o: AlignmentAlgorithm.cpp AlignmentAlgorithm.h IntMatrix.h
	${CC} -c AlignmentAlgorithm.cpp ${OPTS}

CalculatedOvlDistance.o: CalculatedOvlDistance.cpp CalculatedOvlDistance.h
	${CC} -c CalculatedOvlDistance.cpp ${OPTS}

D2.o: D2.cpp D2.h Param.h
	${CC} -c D2.cpp ${OPTS}

ESTAssembly.o: ESTAssembly.cpp ESTAssembly.h Reconstruction.h Graph.h SixTuplesGeneration.h InclusionNodes.h Prim.h
	${CC} -c ESTAssembly.cpp ${OPTS}

Graph.o: Graph.cpp Graph.h Param.h Node.h OvlDistance.h CalculatedOvlDistance.h InclusionNodes.h SixTuple.h Prim.h
	${CC} -c Graph.cpp ${OPTS}

InclusionNodes.o: InclusionNodes.cpp InclusionNodes.h Prim.h
	${CC} -c InclusionNodes.cpp ${OPTS}

OvlDistance.o: OvlDistance.cpp OvlDistance.h Param.h D2.h Alignment.h
	${CC} -c OvlDistance.cpp ${OPTS}

Prim.o: Prim.cpp Prim.h
	${CC} -c Prim.cpp ${OPTS}

Reconstruction.o: Reconstruction.cpp Reconstruction.h SixTuple.h Prim.h Alignment.h InclusionNodes.h Graph.h Param.h SingleBase.h
	${CC} -c Reconstruction.cpp ${OPTS}

SingleBase.o: SingleBase.cpp SingleBase.h
	${CC} -c SingleBase.cpp ${OPTS}

SixTuple.o: SixTuple.cpp SixTuple.h
	${CC} -c SixTuple.cpp ${OPTS}

SixTuplesGeneration.o: SixTuplesGeneration.cpp SixTuplesGeneration.h SixTuple.h Graph.h Param.h InclusionNodes.h
	${CC} -c SixTuplesGeneration.cpp ${OPTS}

clean:
	rm *.o Main