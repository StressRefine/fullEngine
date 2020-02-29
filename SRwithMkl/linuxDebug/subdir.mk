################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../SRPostfF06.cpp \
../SRbasis.cpp \
../SRbasisBrickWedge.cpp \
../SRconstraint.cpp \
../SRcoord.cpp \
../SRedge.cpp \
../SRelemBrickWedge.cpp \
../SRelement.cpp \
../SRerrorCheck.cpp \
../SRface.cpp \
../SRfaceQuad.cpp \
../SRfile.cpp \
../SRforce.cpp \
../SRinput.cpp \
../SRinputUtilities.cpp \
../SRmachDep.cpp \
../SRmap.cpp \
../SRmapBrickWedge.cpp \
../SRmaterial.cpp \
../SRmath.cpp \
../SRmklUtil.cpp \
../SRmodel.cpp \
../SRnode.cpp \
../SRpardiso.cpp \
../SRpardisoSimple.cpp \
../SRpostProcess.cpp \
../SRsaveBreakout.cpp \
../SRsolver.cpp \
../SRstring.cpp \
../SRutil.cpp \
../SRwithMkl.cpp 

OBJS += \
./SRPostfF06.o \
./SRbasis.o \
./SRbasisBrickWedge.o \
./SRconstraint.o \
./SRcoord.o \
./SRedge.o \
./SRelemBrickWedge.o \
./SRelement.o \
./SRerrorCheck.o \
./SRface.o \
./SRfaceQuad.o \
./SRfile.o \
./SRforce.o \
./SRinput.o \
./SRinputUtilities.o \
./SRmachDep.o \
./SRmap.o \
./SRmapBrickWedge.o \
./SRmaterial.o \
./SRmath.o \
./SRmklUtil.o \
./SRmodel.o \
./SRnode.o \
./SRpardiso.o \
./SRpardisoSimple.o \
./SRpostProcess.o \
./SRsaveBreakout.o \
./SRsolver.o \
./SRstring.o \
./SRutil.o \
./SRwithMkl.o 

CPP_DEPS += \
./SRPostfF06.d \
./SRbasis.d \
./SRbasisBrickWedge.d \
./SRconstraint.d \
./SRcoord.d \
./SRedge.d \
./SRelemBrickWedge.d \
./SRelement.d \
./SRerrorCheck.d \
./SRface.d \
./SRfaceQuad.d \
./SRfile.d \
./SRforce.d \
./SRinput.d \
./SRinputUtilities.d \
./SRmachDep.d \
./SRmap.d \
./SRmapBrickWedge.d \
./SRmaterial.d \
./SRmath.d \
./SRmklUtil.d \
./SRmodel.d \
./SRnode.d \
./SRpardiso.d \
./SRpardisoSimple.d \
./SRpostProcess.d \
./SRsaveBreakout.d \
./SRsolver.d \
./SRstring.d \
./SRutil.d \
./SRwithMkl.d 


# Each subdirectory must supply rules for building sources it contributes
%.o: ../%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -I/opt/intel/compilers_and_libraries_2020.0.166/linux/mkl/include -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


