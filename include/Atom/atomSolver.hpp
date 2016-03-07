#ifndef ATOM_SOLVER_HPP
#define ATOM_SOLVER_HPP

#include <vector>
#include <math>
#include <cmath>
#include <stdexcept>

#include <libsgp4/DateTime.h>
#include <libsgp4/Eci.h>
#include <libsgp4/SGP4.h>
#include <libsgp4/Tle.h>

#include "Atom/convertCartesianStateToTwoLineElements.hpp"

namespace atom
{
	template< typename Real, typename Vector3 >
	void executeAtomSolver(
	    const Vector3& departurePosition,
	    const DateTime& departureEpoch,
	    const Vector3& arrivalPosition,
	    const Real timeOfFlight,
	    const Vector3& departureVelocityGuess,
	    Vector3& departureVelocity,
	    Vector3& arrivalVelocity,
	    std::string& solverStatusSummary,
	    int& numberOfIterations,
	    const Tle& referenceTle = Tle( ),
	    const Real earthGravitationalParameter = kMU,
	    const Real earthMeanRadius = kXKMPER,
	    const Real absoluteTolerance = 1.0e-10,
	    const Real relativeTolerance = 1.0e-5,
	    const int maximumIterations = 100 );

	template< typename Real, typename Vector3 >
	void computeNonLinearFunction(
		Vector3 departureVelocity, // input
		Vector3& nonLinearFunctionValue ) // output
	{
		// get target arrival position
		std::vector< Real > targetArrivalPosition( 3 );
		targetArrivalPosition = arrivalPosition;

		// form the departure state vector
		std::vector< Real > departureState( 6 );
		for ( int i = 0; i < 3; i++ )
		{
			departureState[ i ] = departurePosition[ i ];
			departureState[ i + 3 ] = departureVelocity[ i ];
		}

		// convert the departure state vector to TLE
		std::string dummyString = "";
	    int dummyint = 0;
	    const Tle departureTle = convertCartesianStateToTwoLineElements< Real >(
	        departureState,
	        departureEpoch,
	        dummyString,
	        dummyint,
	        referenceTle,
	        earthGravitationalParameter,
	        earthMeanRadius,
	        absoluteTolerance,
	        relativeTolerance,
	        maximumIterations );

	    // propagate departure TLE by time of flight using SGP4 propagator
	    SGP4 sgp4( departureTle );
	    DateTime arrivalEpoch = departureEpoch.AddSeconds( timeOfFlight );
	    Eci arrivalState = sgp4.FindPosition( arrivalEpoch );

	    // compute the value of the non linear function for the corresponding departure velocity
	    nonLinearFunctionValue[ 0 ] = ( arrivalState.Position( ).x - targetArrivalPosition[ 0 ] ) / earthMeanRadius;
	    nonLinearFunctionValue[ 1 ] = ( arrivalState.Position( ).y - targetArrivalPosition[ 1 ] ) / earthMeanRadius;
	    nonLinearFunctionValue[ 2 ] = ( arrivalState.Position( ).z - targetArrivalPosition[ 2 ] ) / earthMeanRadius;
	}

	template< typename Real, typename Vector3, typename Vector2D >
	void computeJacobian(
		Vector3 departureVelocity,
		Vector3 currentFunctionValue,
		Vector2D& jacobian )
	{
		// run the jacobian matrix column loop
		for ( int col = 0; col < 3; col++ )
		{
			Real EPS = 1.0e-8; // approximate square root of machine precision
			Real temp = departureVelocity[ col ];
			Real h = EPS * std::abs( temp );
			if ( h == 0.0 )
			{
				h = EPS;
			} 
			departureVelocity[ col ] = temp + h;
			h = departureVelocity[ col ] - temp;
			std::vector< Real > forwardFunctionValue( 3 );
			computeNonLinearFunction< Real, Vector3 >( departureVelocity, forwardFunctionValue );
			departureVelocity[ col ] = temp; // restore departure velocity value
			// run the jacobian matrix row loop 
			for ( int row = 0; row < 3; row++ )
			{
				jacobian[ row ][ col ] = ( forwardFunctionValue[ col ] - currentFunctionValue[ col ] ) / h;
			}
		}
	}

	// LU decomposition snippet - decomposing the jacobian matrix into lower and upper triangle matrices
	template< typename Real, typename Vector3, typename Vector2D >
	void LUdecomp(
		Vector2D& jacobian, // input and output 
		Vector3& index, // output - records row permutation affected by partial pivoting
		Real& indicator ) // output - indicator whether number of row inter-changes was even or odd 
	{
		const Real tiny = 1.0e-20; // a very small number
		std::vector< Real > rowImplicitScaling( 3 ); // stores implicit scaling of each row of the jacobian matrix
		indicator = 1.0; // no row interchanges yet
		
		// loop over rows to get implicit scaling information for each row of the jacobian matrix.
		// The implicit scaler for each row is the inverse of the largest element in that row.
		for ( int i = 0; i < 3; i++ )
		{
			Real big = 0.0; // stores largest element in the current row 
			for ( int j = 0; j < 3; j++ ) // going over the elements of the current row 'i'
			{
				Real temp = std::abs( jacobian[ i ][ j ] ); // store abs value of a jacobian element
				if ( temp > big )
					big = temp; // looking for a non zero largest element in the current row 'i'
			}
			if ( big == 0.0 ) // no non zero largest element, matrix is singular
				throw std::runtime_error( "ERROR: Singular matrix in routine LUdecomp in atomSolver" );
			rowImplicitScaling[ i ] = 1.0 / big; // store the row scaling, inverse of the largest row element
		}

		int imax = 0;
		int i, j; // indices used in the for loops below
		// loop over the columns of the jacobian matrix (crout's method), lower triangle diagonal elements assumed all unity
		for ( j = 0; j < 3; j++ )
		{
			for ( i = 0; i < j; i++ ) // first procedure, excpt i=j, equation 2.3.12 in Numerical Recipes for c++ 2nd ed. 
			{
				Real sum = jacobian[ i ][ j ]; // storing a[i][j] in 2.3.12
				for ( int k = 0; k < i; k++ )
				{
					sum = sum - jacobian[ i ][ k ] * jacobian[ k ][ j ]; // solving 2.3.12
				}
				jacobian[ i ][ j ] = sum; // storing beta[i][j] in equation 2.3.12
			} // end of first procedure, except for i=j
			
			Real big = 0.0; // to search the largest pivot element (divisor beta) for eqn 2.3.13
			for ( i = j; i < 3; i++ ) // includes i=j of the first procedure for eqn 2.3.12 and i>j of the second procedure for eqn 2.3.13
			{
				Real sum = jacobian[ i ][ j ]; 
				for ( int k = 0; k < j; k++ )
				{
					sum = sum - jacobian[ i ][ k ] * jacobian[ k ][ j ]; // solving RHS bracketed term in equation 2.3.13, pivot division not done
				}
				jacobian[ i ][ j ] = sum;
				Real figureOfMerit = rowImplicitScaling[ i ] * std::abs( sum ); // figure of merit for the pivot element
				if ( figureOfMerit >= big )
				{
					big = figureOfMerit;
					imax = i; // store row value for which the figure of merit for the pivot element is the best
				}
			} // first procedure i=j finished, second procedure i>j partially finished (pivot division still left)
			if ( j != imax ) // for interchanging rows if the condition is true
			{
				for ( int k = 0; k < 3; k++ )
				{
					Real temp = jacobian[ imax ][ k ];
					jacobian[ imax ][ k ] = jacobian[ j ][ k ];
					jacobian[ j ][ k ] = temp;
				}
				indicator = -1.0 * indicator; // change polarity of the row change indicator, indicating a row change happened
				rowImplicitScaling[ imax ] = rowImplicitScaling[ j ]; // interchange row scaling factor as well
			}
			index[ j ] = imax; // stores the row value at which figure of merit for pivot element is the best
			if ( jacobian[ j ][ j ] == 0.0 )
				jacobian[ j ][ j ] = tiny; // changing pivot element from 0.0 to a near zero value, else the jacobian will be singular
			if ( j != 3 - 1 ) // conducting pivot division finally
			{
				Real temp = 1.0 / jacobian[ j ][ j ];
				for ( i = j + 1; i < 3; i++ )
					jacobian[ i ][ j ] *= temp; // division part in equation 2.3.13 finally performed for all i>j, second procedure finished
			} 
		} // loop to next column of the partially modified jacobian matrix
	}

	// LU backward substitution snippet. Solves the multidimensional linear equation A.X = B. 
	// in atom's case, the linear equation is jacobian * dV = -1.0 * nonLinearFunction.
	// foloowing routine adopted from "Numerical recipes in c++" 2nd edition 
	template< typename Real, typename Vector3, typename Vector2D >
	void LUbackSub(
		Vector2D jacobian, // input jacobain matrix in its LU decomposed form
		Vector3 index, // input row permutation vector obtained from the LU decomposition routine
		Vector3& solutionMatrix ) // input as the RHS matrix B in the linear equation and also the output solution matrix X
	{
		int ii = 0; // indexing variable used later for storing the first nonvanishing element of the solution matrix

		// first part is the forward substitution to solve equation 2.3.6
		for ( int i = 0; i < 3; i++ )
		{
			int rowPermutation = index[ i ];
			Real sum = solutionMatrix[ rowPermutation ];
			solutionMatrix[ rowPermutation ] = solutionMatrix[ i ];
			if ( ii != 0 )
			{
				for ( int j = ii - 1; j < i; j++ )
				{
					sum = sum - jacobian[ i ][ j ] * solutionMatrix[ j ];
				}
			}
			else if ( sum != 0.0 )
			{
				ii = i + 1;
			}
			solutionMatrix[ i ] = sum; // store solutions for the forward substitution method
		}

		// second part is the backward substitution method to solve equation 2.3.7
		for ( int i = 3 - 1; i >= 0; i-- )
		{
			Real sum = solutionMatrix[ i ];
			for ( int j = i + 1; j < 3; j++ )
				sum = sum - jacobian[ i ][ j ] * solutionMatrix[ j ];
			solutionMatrix[ i ] = sum / jacobian[ i ][ i ]; // store final solution
		}
	}

	// atom solver snippet
	template< typename Real, typename Vector3 >
	void executeAtomSolver(
	    const Vector3& departurePosition,
	    const DateTime& departureEpoch,
	    const Vector3& arrivalPosition,
	    const Real timeOfFlight,
	    const Vector3& departureVelocityGuess,
	    Vector3& departureVelocity,
	    Vector3& arrivalVelocity,
	    std::string& solverStatusSummary,
	    int& numberOfIterations,
	    const Tle& referenceTle = Tle( ),
	    const Real earthGravitationalParameter,
	    const Real earthMeanRadius,
	    const Real absoluteTolerance,
	    const Real relativeTolerance,
	    const int maximumIterations )
	{
		
	}
} // namespace atom