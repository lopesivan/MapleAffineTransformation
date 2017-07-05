#! /usr/local/bin/maple -q
# -q
(***
 * Test script for the AffineTransformation class
 *)
printlevel := -1;

with (LinearAlgebra);
(* Load some classes *)
read("../modules/affineTransformation.mpl");

# Generates affA
matrixA := Matrix([[1,2,3],[4,5,6]]);
vectorA := Vector([3,3]);

# Generates affB
matrixB := Matrix([[1,2],[-1,2],[8,9]]);
vectorB := Vector([1,-5,2]);

# Input and output vectors for asserts
vectorX := Vector([4,6,9]);
vectorY := Vector([46,103]);
vectorZ := Vector([253, 155, 1297]);

#The lineared form of affA
linearedA := Matrix([[1,2,3,3],[4,5,6,3],[0,0,0,1]]);

#For testing composition
matrixC := matrixA.matrixB;
vectorC := vectorA+(matrixA . vectorB);

#Something square to invert.
matrixD := Matrix([[1,2,3],[1,3,5],[-1,8,9]]);
vectorD := Vector([0,-1,-2]);

matrixE := matrixD;
vectorE := Vector([0,0,0]);

## Build a AffineTransformation with the correct dimensions
try 
	affA := Object ( AffineTransformation, matrixA , vectorA );
	print ("Build Object ... pass.")
catch:
	print ("Build Object ... fail.");
end;

## Get the matrix
try
	if Equal( getLinear( affA ), matrixA ) then
		print ("Get the matrix ... pass.")
	else
		print ("Get the matrix ... fail.")
	fi;
catch: print( "Get the matrix ... fail.")
end;

## Get the vector
try
	if Equal( getTranslation ( affA ), vectorA ) then
		print ("Get the vector ... pass.")
	else
		print ("Get the vector ... fail.")
	fi;
catch:
	print ("Get the vector ... fail.");
end;
## Try to build one with the wrong dimension in the vector
try
	testFailMe := Object ( AffineTransformation, matrixA , Vector([3,3,3]) );
	print ("Prevent wrong dimension ... fail.")
catch:
	print ("Prevent wrong dimension ... pass.");
end;

## Build two AffineTransformations, and multiply them.
try
	affB := Object ( AffineTransformation, matrixB, vectorB );
	affCMaybe := affA.affB;
	#print(AffineTransformation:-getTranslation( affCMaybe ) );
	#print(getLinear(affCMaybe));
	if ( Equal(AffineTransformation:-getTranslation ( affCMaybe ) , vectorC ) 
		and Equal(AffineTransformation:-getLinear ( affCMaybe ), matrixC)
	) then
		print ("Compose two Affine Transformations ... pass.");
	else print ("Compose two Affine Transformations ... fail.");
	fi;
catch:
	print ("Compose two Affine Transformations ... fail.");
end;

## Try to multiply two AffineTransformations of the incorrect dimension
try
	affA.affA;
	print("Prevent Multiplication ... fail.");
catch :
	print("Prevent Multiplication ... pass");
end; #@TODO fix this so that it tests whether or not the error message is correct.

## Get the Dimension of an AffineTransformation, horizontal then vertical.
try
	m,n := AffineDimension( affA );
	if (m = 2 and n = 3) then
		print ("Get dimension of AffineTransformation ... pass.");
	else
		print ("Get dimension of AffineTransformation ... fail.");
	fi
catch:
	print("Get dimension of AffineTransformation ... fail.");
end;

## Get the "Horizontal" Dimension of an AffineTransformation
try
	m := AffineRowDimension( affA );
	if (m = 2) then
		print ("Get row dimension of AffineTransformation ... pass.");
	else
		print ("Get row dimension of AffineTransformation ... fail.");
	fi
catch:
	print ("Get row dimension of AffineTransformation ... fail.");
end;

## Get the "Vertical" Dimension of an AffineTransformation
try
	n := AffineColumnDimension( affA );
	if (n = 3) then
		print ("Get column dimension of AffineTransformation ... pass.");
	else
		print ("Get column dimension of AffineTransformation ... fail.");
	fi
catch:
	print ("Get column dimension of AffineTransformation ... fail.");
end;

## Apply a transformation to a vector.
try
	if Equal( affA.vectorX, vectorY ) then
		print ("Apply transformation to a vector ... pass");
	else
		print ("Apply transformation to a vector ... fail");
	fi;
catch:
	print ("Apply transformation to a vector ... fail");
end;

affD := Object ( AffineTransformation, matrixD, vectorD);
affE := Object (AffineTransformation, matrixE, vectorE);
affF := affD.affE;
affFfromMat := affD.matrixE;
## Multiply an affine transformation and a matrix.
try
	if AffineEqual( affF, affFfromMat ) then
		print("Multiply AffineTransformation and Matrix ... pass");
	else print("Multiply AffineTransformation and Matrix ... fail");
	fi;
catch:
	print("Multiply AffineTransformation and Matrix ... fail");
end;

affF:= affE.affD;
affFfromMat := matrixE.affD;
try
	if AffineEqual( affF, affFfromMat ) then 
		print("Multiply Matrix and AffineTransformation  ... pass");
	else print("Multiply Matrix and AffineTransformation ... fail");
	fi;
catch:
	print("Multiply Matrix and AffineTransformation ... fail");
end;

## Try to apply a transformation to a vector of the incorrect size.
## Try to multiply two AffineTransformations of the incorrect dimension
try
	affA.vectorY;
	print("Prevent Multiplication ... fail.");
catch :
	print("Prevent Multiplication ... pass");
end; #@TODO fix this so that it tests whether or not the error message is correct.

## Evaluate whether or not two Affine Transformations are equal or not.
try 
	affACopy := Object( AffineTransformation, matrixA, vectorA);
	if AffineEqual (affA, affACopy) then
		print ("AffineEqual: Check to see if two affine transformations are equal ... pass.");
	else
		print ("AffineEqual: Check to see if two affine transformations are equal ... fail.");
	fi;
catch:
	print ("AffineEqual: Check to see if two affine transformations are equal ... fail.");
end;

try 
	affACopy := Object( AffineTransformation, matrixA, vectorA);
	if (affA = affACopy) then
		print ("Overloaded =, Check to see if two affine transformations are equal ... pass.");
	else
		print ("Overloaded =, Check to see if two affine transformations are equal ... fail.");
	fi;
catch:
	print ("Overloaded =, Check to see if two affine transformations are equal ... fail.");
end;

try
	if AffineEqual (affA, affB) then
		print ("Check to see if two affine transformations are not equal ... fail.");
	else
		print ("Check to see if two affine transformations are not equal ... pass.");
	fi;
catch:
	print ("Check to see if two affine transformations are not equal ... fail.");
end;

## Get the +1 sized matrix for the affine transformation
#try
	# http://negativeprobability.blogspot.com/2011/11/affine-transformations-and-their.html
#	if Equal(MakeLinear( affA ),linearedA) then
#		print("MakeLinear ... pass.");
#	else
#		print("MakeLinear ... fail.");
#	fi;
#catch :
#	print("MakeLinear ... fail.")
#end;

## Build an affine transformation from a single, "+1" matrix.
#try
#	affACopy := CreateFromLargerLinear( AffineTransformation, linearedA );
#	if Equal(linearedA, MakeLinear(affACopy)) then
#		print ( "CreateFromLargerLinear ... pass." );
#	else
#		print ( "CreateFromLargerLinear ... fail.");
#	fi;
#catch :
#	print("CreateFromLargerLinear ... fail.");
#end;

## Try to build an affine transformation from an "incorrect", "+1" matrix.
#   not sure if this is right...or what we want to do.  Is there something
#	"like" the Jordan-normal form, where the lower-right entry is a one,
#	and all other entries on the bottom row are zero?

## Check to see if a matrix is an affine transformation.
#try
#	if isLargeAffine ( AffineTransformation, linearedA ) then
#		print ( "isLargeAffine ... pass." );
#	else
#		print ( "isLargeAffine ... fail." );
#	fi;
#catch :
#	print("isLargeAffine ... fail.");
#end;

## Invert an affine transformation
try
	affD := Object ( AffineTransformation, matrixD, vectorD);
	matInvD := MatrixInverse(matrixD);
	vectInvD := (-1 * matInvD) . vectorD; 
	#whattype(invD); print(whattype(affD));
	affDinv := affD:-AffineInverse ( affD );
	#affDinv:-getLinear(affDinv);  affDinv:-getTranslation(affDinv);
	if (
		Equal(affDinv:-getLinear (affDinv), matInvD ) 
		and Equal(affDinv:-getTranslation(affDinv), vectInvD)
	)
	then print("AffineInverse ... pass.");
	else print("AffineInverse ... fail.");
	fi;
catch :
	print("AffineInverse ... fail.");
end; #@TODO calculate the inverse of that thing, and check to see if it worked.

try
	affD:= Object(AffineTransformation,matrixD,vectorD);
	affE:= affD.(affD.affD);
	affEtest := affD:-AffineSelfMultiply(affD,3);
	if (affE = affEtest) then
		print ("AffineSelfMultiply ... pass.");
	else print ("AffineSelfMultiply ... fail.");
	fi;
catch :
	print ("AffineSelfMultiply ... fail.");
end;

## Try the overloaded carrot operation
try 
	affD:= Object(AffineTransformation,matrixD,vectorD);
	affE:= affD.(affD.affD);
	affEtest := affD^3;
	if (affE = affEtest) then
		print ("Carrot positive ... pass.");
	else print ("Carrot positive ... fail.");
	fi;
catch :
	print ("Carrot positive ... fail.");
end;	

try 
	affD:= Object(AffineTransformation,matrixD,vectorD);
	affDinv := affD:-AffineInverse ( affD );
	affE := affDinv . affDinv;
	affEtest := affD^(-2);
	if (affE = affEtest) then
		print ("Carrot negative ... pass.");
	else print ("Carrot negative ... fail.");
	fi;
catch :
	print ("Carrot negative ... fail.");
end;

try
	affD:= Object(AffineTransformation,matrixD,vectorD);
	affE := Object ( AffineTransformation, Matrix([[1,0,0],[0,1,0],[0,0,1]]),Vector([0,0,0]));
	affEtest := affD^(0);
	if (affE = affEtest) then
		print ("Carrot zero ... pass.");
	else print ("Carrot zero ... fail.");
	fi;
catch :
	print ("Carrot zero ... fail.");
end;

# Create an affine transformation with vector, matrix.
try
	affAtest:= Object (AffineTransformation, vectorA, matrixA);
	if (affAtest = affA) then
		print ("Construct in reverse order ... pass");
	else ("Construct in reverse order ... fail");
	fi;
catch :
	print ("Construct in reverse order ... fail");
end;
# Create an affine transformation with overloaded function caller
try
	affAtest := AffineTransformation( matrixA, vectorA );
	if (affAtest = affA) then
		print ("Call as function ... pass");
	else ("Call as function ... pass");
	fi;
catch :
	print ("Call as function ... pass");
end;	