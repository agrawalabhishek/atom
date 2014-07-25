/*    
 * Copyright (c) 2010-2014, Delft University of Technology
 * Copyright (c) 2010-2014, K. Kumar (me@kartikkumar.com)
 * All rights reserved.
 * See http://bit.ly/12SHPLR for license details.
 */

#include <vector>

#include <Eigen/Core>
 
#include <gsl/gsl_vector.h>

#include <libsgp4/Tle.h>

namespace atom
{

//! Typedef for Vector6d.
typedef Eigen::Matrix< double, 6, 1 > Vector6d;

//! Convert Cartesian state to TLE (Two Line Elements).
/*!
 * Converts a given Cartesian state (position, velocity) to an equivalent TLE. 
 *
 * This function makes use of a root-finder to solve a non-linear system. Locating the root of the 
 * non-linear system corresponds finding the TLE that, when evaluated at its epoch using the 
 * SGP4/SDP4 propagator, yields the target Cartesian state (within tolerance).
 *
 * Details of the underlying non-linear system and algorithm are catalogued by Kumar, et al. (2014).
 *
 * \param cartesianState Cartesian state.
 * \param epoch Epoch associated with Cartesian state.
 * \param earthGravitationalParameter Earth gravitational parameter [m^3 s^-2].
 * \param referenceTle Reference Two Line Elements. 
 * \param tolerance Tolerance used to check if root-finder has converged [default: 1.0e-8].
 * \param isPrintProgress Should the progress of the root-finder be printed to console? This is 
            useful for debugging [default: false].
 * \return TLE object that generates target Cartesian state when propagated with SGP4 propagator to
 *           target epoch.
 * \sa evaluateCartesianToTwoLineElementsSystem()
 */
const Tle convertCartesianStateToTwoLineElements( const Vector6d cartesianState,
                                                  const DateTime epoch,
                                                  const double earthGravitationalParameter,
                                                  const Tle referenceTle,                                                  
                                                  const double tolerance = 1.0e-8,
                                                  const bool isPrintProgress = false );

//! Evaluate system of non-linear equations for converting Cartesian state to TLE.
/*!
 * Evaluates system of non-linear equations to find TLE corresponding with target Cartesian 
 * state. The system of non-linear equations used is:
 *  \f[ 
 *      F = 0 = \bar{x}^{new} - \bar{x}^{target}
 *  \f]
 * where \f$\bar{x}^{new}\f$ is the new Cartesian state computed by updating the TLE mean elements 
 * and propagating the TLE using the SGP4 propagator, and \f$\bar{x}^{target}\f$ is the target
 * Cartesian state.
 * \param indepedentVariables Vector of independent variables used by the root-finder.
 * \param parameters Parameters required to compute the objective function.
 * \param functionValues Vector of computed function values. 
 */
int evaluateCartesianToTwoLineElementsSystem( const gsl_vector* independentVariables, 
                                              void* parameters, 
                                              gsl_vector* functionValues );

//! Update TLE mean elements.
/*!
 * Updates mean elements stored in TLE based on current osculating elements. This function 
 * uses the osculating elements to replace mean elements (converts units and computes mean anomaly
 * and mean motion).
 * \param newKeplerianElements New Keplerian elements.
 * \param oldTle TLE in which the mean elements are to be replaced.
 * \param earthGravitationalParameter Earth gravitational parameter [m^3 s^-2].
 * \return New TLE with mean elements updated.
 */
const Tle updateTleMeanElements( const Eigen::VectorXd newKeplerianElements, 
    							 const Tle oldTle,
    							 const double earthGravitationalParameter );

//! Container of parameters used by Cartesian-to-TLE objective function.
struct CartesianToTwoLineElementsObjectiveParameters
{ 
public:

    //! Constructor taking parameter values.
    CartesianToTwoLineElementsObjectiveParameters( 
        const Vector6d aTargetState,
        const double anEarthGravitationalParameter,
        const Tle someReferenceTle )
        : targetState( aTargetState ),
          earthGravitationalParameter( anEarthGravitationalParameter ),
 		  referenceTle( someReferenceTle )
    { }

    //! Target state in Cartesian elements.
    const Eigen::VectorXd targetState;

    //! Earth gravitational parameter [kg m^3 s^-2].
    const double earthGravitationalParameter;

    //! Reference TLE.
    const Tle referenceTle;

protected:

private:
};

} // namespace atom