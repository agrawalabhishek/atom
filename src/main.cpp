#include <iostream>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

#include <libsgp4/DateTime.h>
#include <libsgp4/Eci.h>
#include <libsgp4/Globals.h>
#include <libsgp4/SGP4.h>
#include <libsgp4/Tle.h>
#include <libsgp4/Vector.h>

#include <Astro/astro.hpp>
#include <SML/sml.hpp>

#include "Atom/convertCartesianStateToTwoLineElements.hpp"

typedef std::vector< double > Vec;

int main()
{
    // std::cout << "Hello!" << std::endl;
    // std::cout << std::endl;
    std::cout.precision( 15 );

    //
    const double apsisMinimum = 6700.0;
    const double apsisMaximum = 50000.0;

    // const double inclinationMinimum = 0.0;
    // const double inclinationMaximum = 180.0;

    // const double argumentOfPeriapsisMinimum = 0.0;
    // const double argumentOfPeriapsisMaximum = 360.0;

    // const double longitudeOfAscendingNodeMinimum = 0.0;
    // const double longitudeOfAscendingNodeMaximum = 360.0;

    // const double trueAnomalyMinimum = 0.0;
    // const double trueAnomalyMaximum = 360.0;

    const int trials = 100000;

    //
    std::random_device rand_seed;
    std::default_random_engine random_generator( rand_seed( ) );

    std::uniform_real_distribution< double > uniformApsis( apsisMinimum, apsisMaximum );
    std::uniform_real_distribution< double > uniformAngles( 0.0, 360.0 );

    DateTime epoch( 2014, 6, 13, 14, 59, 21 );
    // std::cout << epoch << std::endl;
    // std::cout << std::endl;

    int countUndefinedExceptions = 0;
    int countDecayedExceptions = 0;
    int countPlExceptions = 0;
    int countElsqExceptions = 0;
    int countIterationsExceptions = 0;
    int countFindPositionExceptions = 0;
    int countXPositionRelativeErrorExceptions = 0;
    int countYPositionRelativeErrorExceptions = 0;
    int countZPositionRelativeErrorExceptions = 0;
    int numberOfIterations = 0;
    std::string solverStatusSummary;

    for ( unsigned int trial = 0; trial < trials; ++trial )
    {
        Vec keplerianElements( 6 );
        const double apsis1 = uniformApsis( random_generator );
        const double apsis2 = uniformApsis( random_generator );

        // std::cout << apsis1 << ", " << apsis2 << std::endl;

        double periapsis = apsis1;
        if ( apsis1 > apsis2 )
        {
            periapsis = apsis2;
        }

        keplerianElements[ astro::eccentricityIndex ] = std::fabs( apsis1 - apsis2 ) / ( apsis1 + apsis2 );
        // keplerianElements[ astro::eccentricityIndex ] = 0.0001;
        keplerianElements[ astro::semiMajorAxisIndex ] = periapsis / ( 1.0 - keplerianElements[ astro::eccentricityIndex ] );
        keplerianElements[ astro::inclinationIndex ] = sml::convertDegreesToRadians( sml::computeModulo( uniformAngles( random_generator ), 180.0 ) );
        keplerianElements[ astro::argumentOfPeriapsisIndex ] = sml::convertDegreesToRadians( uniformAngles( random_generator ) );
        keplerianElements[ astro::longitudeOfAscendingNodeIndex ] = sml::convertDegreesToRadians( uniformAngles( random_generator ) );
        keplerianElements[ astro::trueAnomalyIndex ] = sml::convertDegreesToRadians( uniformAngles( random_generator ) );

        // for ( unsigned int i = 0; i < 6; ++i )
        // {
        //     std::cout << keplerianElements[ i ] << std::endl;
        // }
        // std::cout << std::endl;

        Vec cartesianState = astro::convertKeplerianToCartesianElements( keplerianElements, kMU );
        // for ( unsigned int i = 0; i < 6; ++i )
        // {
        //     std::cout << cartesianState[ i ] << std::endl;
        // }
        // std::cout << std::endl;

        Tle virtualTle;
        try
        {
            virtualTle = atom::convertCartesianStateToTwoLineElements< double >( cartesianState, epoch, solverStatusSummary, numberOfIterations );
        }
        catch( std::exception& e )
        {
            // std::cout << e.what( ) << std::endl;
            std::cout << solverStatusSummary << std::endl;

            if ( strcmp( e.what( ), "Error: Satellite decayed" ) == 0 )
            {
                ++countDecayedExceptions;
                // std::cout << e.Position( ) << std::endl;
                // std::cout << e.Velocity( ) << std::endl;
                // std::cout << e.Decayed( ) << std::endl;
            }

            else if ( strcmp( e.what( ), "Error: (pl < 0.0)" ) == 0 )
            {
                ++countPlExceptions;
            }

            else if ( strcmp( e.what( ), "Error: (elsq >= 1.0)" ) == 0 )
            {
                ++countElsqExceptions;
            }

            else
            {
                ++countUndefinedExceptions;
            }

            // for ( unsigned int i = 0; i < 5; ++i )
            // {
            //     std::cout << keplerianElements[ i ] << ", ";
            // }
            // std::cout << keplerianElements[ 5 ] << std::endl;
            // std::cout << std::endl;
            continue;
        }

        if ( numberOfIterations > 100 )
        {
            throw( "iterations > 100 " );
            ++countIterationsExceptions;
            continue;
        }

        // std::cout << virtualTle << std::endl;
        // std::cout << std::endl;

        Vector position( 0.0, 0.0, 0.0, 0.0 );
        Vector velocity( 0.0, 0.0, 0.0, 0.0 );
        Eci state( epoch, position, velocity );
        try
        {
            SGP4 sgp4( virtualTle );
            state = sgp4.FindPosition( 0.0 );
        }
        catch( std::exception& e )
        {
            ++countFindPositionExceptions;
            std::cout << e.what( ) << std::endl;
        }

        const double tolerance = 1.0e-5;

        const double xPositionRelativeError = ( state.Position( ).x - cartesianState[ 0 ] ) / cartesianState[ 0 ];
        if ( xPositionRelativeError > tolerance )
        {
            std::cout << "x-position relative error: "
                      << xPositionRelativeError << ", "
                      << state.Position( ).x  << ", "
                      << cartesianState[ 0 ] << std::endl;
            // for ( unsigned int i = 0; i < 5; ++i )
            // {
            //     std::cout << keplerianElements[ i ] << ", ";
            // }
            // std::cout << keplerianElements[ 5 ] << std::endl;
            // std::cout << std::endl;
            ++countXPositionRelativeErrorExceptions;
        }

        const double yPositionRelativeError = ( state.Position( ).y - cartesianState[ 1 ] ) / cartesianState[ 1 ];
        if ( yPositionRelativeError > tolerance )
        {
            std::cout << "y-position relative error: "
                      << yPositionRelativeError << ", "
                      << state.Position( ).y << ", "
                      << cartesianState[ 1 ]<< std::endl;
            // for ( unsigned int i = 0; i < 5; ++i )
            // {
            //     std::cout << keplerianElements[ i ] << ", ";
            // }
            // std::cout << keplerianElements[ 5 ] << std::endl;
            // std::cout << std::endl;
            ++countYPositionRelativeErrorExceptions;
        }

        const double zPositionRelativeError = ( state.Position( ).z - cartesianState[ 2 ] ) / cartesianState[ 2 ];
        if ( zPositionRelativeError > tolerance )
        {
            std::cout << "z-position relative error: "
                      << zPositionRelativeError << ", "
                      << state.Position( ).z << ", "
                      << cartesianState[ 2 ] << std::endl;
            // for ( unsigned int i = 0; i < 5; ++i )
            // {
            //     std::cout << keplerianElements[ i ] << ", ";
            // }
            // std::cout << keplerianElements[ 5 ] << std::endl;
            // std::cout << std::endl;
            ++countZPositionRelativeErrorExceptions;
        }
    }

    const int countExceptions = countDecayedExceptions
                              + countPlExceptions
                              + countElsqExceptions
                              + countUndefinedExceptions
                              + countIterationsExceptions
                              + countFindPositionExceptions
                              + countXPositionRelativeErrorExceptions
                              + countYPositionRelativeErrorExceptions
                              + countZPositionRelativeErrorExceptions;

    std::cout << std::endl;
    std::cout << "Number of Decayed exception cases: " << countDecayedExceptions << std::endl;
    std::cout << "Number of Pl exception cases: " << countPlExceptions << std::endl;
    std::cout << "Number of Elsq exception cases: " << countElsqExceptions << std::endl;
    std::cout << "Number of undefined exception cases: " << countUndefinedExceptions << std::endl;
    std::cout << "Number of iterations exception cases: " << countIterationsExceptions << std::endl;
    std::cout << "Number of FindPosition exception cases: " << countFindPositionExceptions << std::endl;
    std::cout << "Number of x-position relative error exception cases: " << countXPositionRelativeErrorExceptions << std::endl;
    std::cout << "Number of y-position relative error exception cases: " << countYPositionRelativeErrorExceptions << std::endl;
    std::cout << "Number of z-position relative error exception cases: " << countZPositionRelativeErrorExceptions << std::endl;
    std::cout << "Number of exception cases: " << countExceptions << std::endl;
    std::cout << std::endl;

    return 0;
}