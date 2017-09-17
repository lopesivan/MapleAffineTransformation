(***** AffineTransformation Object
 * @author Scott Grizzard <sgrizzar@mail.usf.edu>
 *
 * An affine transformation A from an m-dimensional space X to an n-dimensional space Y
 * is an matrix M (with horizontal dimension m and vertical dimension n), and an n-vector v,
 * so that A(x) = Mx + v.
 * 
 * Composition (multiplication) of two affine transformations is read
 * right to left (as matrix multiplication is).
 *
**)

module AffineTransformation() option object;

	description "An affine transformation A from an m-dimensional space X to an n-dimensional space Y is an matrix M (with horizontal dimension m and vertical dimension n), and an n-vector v, so that A(x) = Mx + v.";
	
	local
		linearPortion::Matrix,		#Matrix holding the linear portion of the Affine Transformation
		translationVector::Vector;  #Vector holding the translation portion of the Affine Transformation

	export 
		(**
	 	 * Variable gets.
	 	 *)	
		getLinear::static := proc ( self::AffineTransformation ) :: Matrix;
			description "Returns the linear portion of the transformation in matrix form";
			self:-linearPortion;
		end proc,
		getTranslation::static := proc (self::AffineTransformation ) :: Vector;
			description "Returns the translation vector portion of the transformation in vector form";
			self:-translationVector;
		end proc,
		(**
		 * Constructor method	
		 *)
		ModuleCopy::static := overload ([
			proc ( self::AffineTransformation,
				proto::AffineTransformation, 
				linearPortion::Matrix,
				translationVector::Vector) :: AffineTransformation;
				
				description "Constructor for AffineTransformation object.  Takes the linear portion as a matrix, and the translation vector as a vector.";
				option overload;
				
				#@TODO add checkers to make sure these are of the proper size.
				if self:-dimensionChecker(linearPortion, translationVector) then
					self:-linearPortion := linearPortion;
					self:-translationVector := translationVector;
				else error ("Row dimension of matrix portion of transformation does not agree with dimension of translation vector.");
				fi;
			end proc,
			
			proc ( self::AffineTransformation,
				proto::AffineTransformation,
				translationVector::Vector, 
				linearPortion::Matrix) :: AffineTransformation;
				
				description "Constructor for AffineTransformation object.  Takes the linear portion as a matrix, and the translation vector as a vector.";
				option overload;
				
				#@TODO add checkers to make sure these are of the proper size.
				if self:-dimensionChecker(linearPortion, translationVector) then
					self:-linearPortion := linearPortion;
					self:-translationVector := translationVector;
				else error ("Row dimension of matrix portion of transformation does not agree with dimension of translation vector.");
				fi;
			end proc
			]),
		(**
		 * Apply method ... if the object is called as a function
		 *)
		ModuleApply::static := proc( )
        	Object( AffineTransformation , _passed );
    	end proc,
    
		(**
		 * Pretty printing
		 *)
		ModulePrint::static := proc (self::AffineTransformation) :: list;
			[self:-translationVector,self:-linearPortion];
		end proc,
		
		(**
		 * Simplify overloaded
		 * @TODO Make this (or something like it) work!
		 *)
		(**simplify::static := overload([
			proc(self::AffineTransformation, $)::AffineTransformation;
				option overload;
				self:-linearPortion := simplify(self:-linearPortion);
				self:-translationVector := simplify(self:-translationVector);
				return (self);
			end proc
		]),*)
		(**
		 * Composes two affine transformations.
		 *
		 *)
		AffineAffineMultiply::static := proc ( affA::AffineTransformation, affB::AffineTransformation, $ ) :: AffineTransformation;
			description "Composes two affine transformations, in the same order as matrix multiplication (reading right to left).";
			return Object ( AffineTransformation , 
				LinearAlgebra[MatrixMatrixMultiply](getLinear(affA),getLinear(affB)),
				getTranslation(affA) + LinearAlgebra[MatrixVectorMultiply] (getLinear(affA),getTranslation(affB) )
			);
		end proc,
		(**
		 * Applies an affine transformation to a vector
		 *)
		AffineVectorMultiply::static := proc ( affA::AffineTransformation, vectX::Vector, $ )
			::Vector;
			simplify(LinearAlgebra[MatrixVectorMultiply] ( getLinear( affA ) , vectX ) + getTranslation (affA));
		end proc,
		(**
		 * Lets you multiply a matrix and an affine transformation.
		 *)
		AffineMatrixMultiply::static := proc (affA::AffineTransformation, matM::Matrix, $ )
			:: AffineTransformation;
			return Object ( AffineTransformation ,
				simplify(LinearAlgebra[MatrixMatrixMultiply](getLinear(affA), matM )),
				getTranslation(affA)  #Plus zero vector
			);
		end proc,
		MatrixAffineMultiply::static := proc (matM::Matrix, affA::AffineTransformation, $ )
			:: AffineTransformation;
			return Object ( AffineTransformation,
				LinearAlgebra[MatrixMatrixMultiply](matM, getLinear(affA)),
				LinearAlgebra[MatrixVectorMultiply](matM, getTranslation(affA))
			);
		end proc,
		ScalerAffineMultiply::static := proc (c::numeric, affA::AffineTransformation, $ )
			:: AffineTransformation;
			local matC::Matrix;
			matC := LinearAlgebra[MatrixScalarMultiply](LinearAlgebra[IdentityMatrix](affA:-AffineRowDimension(affA),affA:-AffineColumnDimension(affA)),c);
			return (affA:-MatrixAffineMultiply(matC,affA));
		end proc,
		
		AffineScalerMultiply::static := proc (affA::AffineTransformation, c::numeric, $)
			:: AffineTransformation;
			local matC::Matrix;
			matC := LinearAlgebra[MatrixScalarMultiply](LinearAlgebra[IdentityMatrix](affA:-AffineRowDimension(affA),affA:-AffineColumnDimension(affA)),c);
			return (affA:-AffineMatrixMultiply(affA,matC));
		end proc,
		
		AffineSelfMultiply::static := proc (affA::AffineTransformation, power::posint) :: AffineTransformation;
			description "Multiplies an affine transformation by itself power times.";
			local newAffine::AffineTransformation,
				newPower::posint;
			
			newAffine := affA;
			newPower:=power;
			
			while (newPower > 1) do
				newAffine:= affA:-AffineAffineMultiply(newAffine, affA);
				newPower := newPower - 1;
			od;
			return newAffine;
		end proc,
		(**
		 * Overloaded composition operator
		 *)
		`.`::static := overload(
			[	
				proc ( affA::AffineTransformation, affB::AffineTransformation, $ ) :: AffineTransformation;
					option overload;
					return AffineAffineMultiply( affA, affB);
				end proc,
				
				proc ( affA::AffineTransformation, vectX::Vector, $ ) :: Vector;
					option overload;
					return AffineVectorMultiply (affA, vectX);
				end proc,
				
				proc ( affA::AffineTransformation, matM::Matrix, $ ) :: AffineTransformation;
					option overload;
					return affA:-AffineMatrixMultiply (affA, matM);
				end proc,
				
				proc ( matM::Matrix, affA::AffineTransformation, $) :: AffineTransformation;
					option overload;
					return affA:-MatrixAffineMultiply (matM, affA);
				end proc,
				
				proc ( c::numeric, affA::AffineTransformation ) :: AffineTransformation;
					option overload;
					return affA:-ScalerAffineMultiply(c,affA);
				end proc,
				
				proc (affA::AffineTransformation, c::numeric, $) :: AffineTransformation;
					option overload;
					return affA:-AffineScalerMultiply(affA,c);
				end proc
				#,
				#proc (a,b,$)
				#	option overload;
				#	error("got bad stuff", type(a, AffineTransformation), type(b, AffineTransformation));
				#end proc
			]
		),
		(**
		 * Inverts an affine transformation, returning a new AffineTransformation.
		 *)
		AffineInverse::static := proc ( self::AffineTransformation ) :: AffineTransformation;
			local matInv :=  simplify(LinearAlgebra[MatrixInverse]( self:-linearPortion ));
			#print (matInv);
			local neg_matInv := LinearAlgebra[MatrixScalarMultiply](matInv,-1);
			#print (neg_matInv);
			local vectInv := LinearAlgebra[MatrixVectorMultiply](neg_matInv,self:-translationVector);
			#print (vectInv);
			return Object ( AffineTransformation,
				matInv,
				vectInv
			);
		end proc,
		#@TODO overload carrot to do powers and inverses.
		(**
		 * Dimension functions
		 *)
		AffineDimension::static := proc ( self::AffineTransformation ) #what do you call the thing this returns?
			description "Returns the dimension of the matrix containing the linear portion of the transformation";
			return LinearAlgebra[Dimension](self:-linearPortion);
		end proc,
		AffineRowDimension::static := proc (self::AffineTransformation ) :: int;
			description "Returns the row dimension of the affine transformation.";
			return LinearAlgebra[RowDimension](self:-linearPortion);
		end proc,
		AffineColumnDimension::static := proc ( self::AffineTransformation ) :: int;
			description "Returns the column dimension of the affine transformation.";
			return LinearAlgebra[ColumnDimension](self:-linearPortion);
		end proc,
		(**
		 * Equal Function
		 *)
		AffineEqual::static := proc ( affA::AffineTransformation, affB::AffineTransformation ) :: bool;
		 	description "Checks to see if two Affine Transformations are equal.";
		 	if (
		 		LinearAlgebra[Equal](affA:-linearPortion,affB:-linearPortion) and 
		 		LinearAlgebra[Equal](affA:-translationVector, affB:-translationVector)
		 	) then
		 		return true;
		 	else
		 		return false;
		 	fi;
		end proc,
		`=`::static := overload([
			proc (affA::AffineTransformation, affB::AffineTransformation, $) :: bool;
				option override;
				return AffineTransformation:-AffineEqual (affA, affB);
			end proc
		]),
		(**
		 *	Overloading the power operator.
		 *  If you get a one, then just send back the affine transformation.  Otherwise, make
		 *	sure the transformation is square.  If so, then go through the following switch:
		 *		power is positive, in which case, apply AffineSelfMultiply(affA, power).
		 *		Power is negative, in which case, invert affA, and send the inverse to AffineSelfMultiply,
		 *			with -power.
		 *		Power is zero, in which case, return the identity transformation (identity matrix and zero vector).
		 *)
		`^`::static := overload([
			proc (affA::AffineTransformation, power::integer, $) :: AffineTransformation ;
				option overload;
				#check if one.  if so, return affA.
				if (power = 1) then return affA; fi;
				
				#check square.  if not, raise exception
				if (affA:-checkSquare(affA:-getLinear(affA))) then
					# series of ifs
					if ( power = 0 ) then
						#return (Object (AffineTransformation, 
						#	LinearAlgebra[IndentityMatrix](affA:-AffineRowDimension(affA),affA:-AffineColumnDimension(affA)),
						#	LinearAlgebra[ZeroVector](affA:-AffineRowDimension(affA))	
						#	));
						#error("ugly");
						#Case is never read because of the way Maple works...it will always the scalar 1
						elif ( power > 0 ) then
						return (affA:-AffineSelfMultiply(affA,power));
					elif (power < 0 ) then
						return (affA:-AffineSelfMultiply(affA:-AffineInverse(affA), -1 * power));
					fi;
				else error("Non-square Affine Transformation" );
				fi;
			end proc
		]);
		# @TODO Overload the ^ operator to take integers and give inverses or powers.
				
	
	local
		(**
		 * Check to make sure the vector dimension matches the dimension of the matrix.
		 *)
		dimensionChecker::static := proc ( linearPortion::Matrix, translationVector::Vector ) :: bool;
			description "Checks whether or not the row dimension of the matrix is the same as the dimension of the translation vector.";
			if LinearAlgebra[RowDimension](linearPortion) = LinearAlgebra[Dimension](translationVector) then
				return true;
			else return false;
			fi;
		end proc,
		
		(**
		 * Check to make sure you have a square matrix.
		 *)
		checkSquare::static := proc ( linearPortion::Matrix ) :: bool;
			return (
				LinearAlgebra[RowDimension](linearPortion) = LinearAlgebra[ColumnDimension](linearPortion)
			);
		end proc;

end module;