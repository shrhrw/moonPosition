/******************************************************************************
 * Copyright (c) 2021 ispace technologies US
 *
 * Reference(s):
 * (1) Meeus,Jean."Astronomical Algorithms",
 * 2nd Edition, Willmann-Bell, Inc.
 * Chapter 22. (Nutation and the Obliquity of the Ecliptic)
 * Chapter 47. (Position of the Moon)
 ******************************************************************************
*/

#include "../include/moonPosition.h"

constexpr double MoonPosition::LONGITUDE_PERIODIC_TERM_MATRIX[60][4];
constexpr double MoonPosition::LATITUDE_PERIODIC_TERM_MATRIX[60][4];
constexpr double MoonPosition::LONGITUDE_SINE_COEFF[60];
constexpr double MoonPosition::LONGITUDE_COSINE_COEFF[60];
constexpr double MoonPosition::LATITUDE_SINE_COEFF[60];

MoonPosition::MoonPosition() :
    longitudeIntegrated(),
    latitudeIntegrated(),
    earthToMoonRadius(),
    dynamicalTime(),
    earthEccentricity(),
    actionOfVenus(),
    actionOfJupiter(),
    actionOfSomething(),
    moonMeanElongation(),
    moonMeanLongitude(),
    moonArgumentOfLatitude(),
    moonMeanAnomaly(),
    sunMeanAnomaly(),
    geocentricLongitude(),
    geocentricLatitude(),
    earthToMoonDistance(),
    julianEphemerisDay(),
    cartesianCoordinates()
{}

MoonPosition::~MoonPosition()
{}

double MoonPosition::Reduce( double value )
{
    value = fmod( value, 360.0 );
    return ( value < 0 ? value + 360.0 : value );     
}

// Calculate Dynamical Time (T) from Julian Ephemeris Day (JDE).
// JDE is the number of centuries elapsed since the Julian Epoch (J2000).
void MoonPosition::CalcDynamicTime()
{
    this->dynamicalTime = ( julianEphemerisDay - DYNAMIC_TIME_COEFFS[0] ) /
                        DYNAMIC_TIME_COEFFS[1];
}

// Moon mean longitude, referred to as the mean equinox of the date
// including constant term of the effect of light-time (-0``.70) (L`)(degrees)
void MoonPosition::CalcMoonMeanLongitude()
{
    this->moonMeanLongitude = Reduce( ( MOON_MEAN_LONGITUDE_COEFFS[0] +
                                      MOON_MEAN_LONGITUDE_COEFFS[1] * dynamicalTime -
                                      MOON_MEAN_LONGITUDE_COEFFS[2] * pow( dynamicalTime, 2 ) +
                                      MOON_MEAN_LONGITUDE_COEFFS[3] * pow( dynamicalTime, 3 ) -
                                      MOON_MEAN_LONGITUDE_COEFFS[4] * pow( dynamicalTime, 4 ) ) );
}

// Mean elongation of the Moon (D)(degrees)
void MoonPosition::CalcMoonMeanElongation()
{
    this->moonMeanElongation = Reduce( ( MOON_MEAN_ELONGATION_COEFFS[0] +
                                       MOON_MEAN_ELONGATION_COEFFS[1] * dynamicalTime -
                                       MOON_MEAN_ELONGATION_COEFFS[2] * pow( dynamicalTime, 2 ) +
                                       MOON_MEAN_ELONGATION_COEFFS[3] * pow( dynamicalTime, 3 ) -
                                       MOON_MEAN_ELONGATION_COEFFS[4] * pow( dynamicalTime, 4 ) ) );
}

// Sun mean anomaly (M)(degrees)
void MoonPosition::CalcSunMeanAnomaly()
{
    this->sunMeanAnomaly = Reduce( ( SUN_MEAN_ANOMALY_COEFFS[0] +
                                   SUN_MEAN_ANOMALY_COEFFS[1] * dynamicalTime -
                                   SUN_MEAN_ANOMALY_COEFFS[2] * pow( dynamicalTime, 2 ) +
                                   SUN_MEAN_ANOMALY_COEFFS[3] * pow( dynamicalTime, 3 ) ) );
}

// Moon mean anomaly (M`)(degrees)
void MoonPosition::CalcMoonMeanAnomaly()
{
    this->moonMeanAnomaly = Reduce( ( MOON_MEAN_ANOMALY_COEFFS[0] +
                                    MOON_MEAN_ANOMALY_COEFFS[1] * dynamicalTime +
                                    MOON_MEAN_ANOMALY_COEFFS[2] * pow( dynamicalTime, 2 ) +
                                    MOON_MEAN_ANOMALY_COEFFS[3] * pow( dynamicalTime, 3 ) -
                                    MOON_MEAN_ANOMALY_COEFFS[4] * pow( dynamicalTime, 4 ) ) );
}

// Moon argument of latitude (mean distance of the Moon from its ascending node) (F)(degrees)
void MoonPosition::CalcMoonArgumentOfLatitude()
{
    this->moonArgumentOfLatitude = Reduce( ( MOON_ARG_OF_LATITUDE_COEFFS[0] +
                                           MOON_ARG_OF_LATITUDE_COEFFS[1] * dynamicalTime -
                                           MOON_ARG_OF_LATITUDE_COEFFS[2] * pow( dynamicalTime, 2 ) -
                                           MOON_ARG_OF_LATITUDE_COEFFS[3] * pow( dynamicalTime, 3 ) +
                                          MOON_ARG_OF_LATITUDE_COEFFS[4] * pow( dynamicalTime, 4 ) ) );
}

void MoonPosition::CalcActionOfVenus()
{
    this->actionOfVenus = Reduce( ( ACTION_OF_VENUS_COEFFS[0] + ACTION_OF_VENUS_COEFFS[1] * dynamicalTime ) );
}

void MoonPosition::CalcActionOfJupiter()
{
    this->actionOfJupiter = Reduce( ( ACTION_OF_JUPITER_COEFFS[0] + ACTION_OF_JUPITER_COEFFS[1] * dynamicalTime ) );
}

void MoonPosition::CalcActionOfSomething()
{
    this->actionOfSomething = Reduce( ( ACTION_OF_SOMETHING_COEFFS[0] + ACTION_OF_SOMETHING_COEFFS[1] * dynamicalTime ) );
}

void MoonPosition::CalcEarthEccentricity()
{
    this->earthEccentricity = EARTH_ECCENTRICITY_COEFFS[0] -
                              EARTH_ECCENTRICITY_COEFFS[1] * dynamicalTime -
                              EARTH_ECCENTRICITY_COEFFS[2] * pow( dynamicalTime, 2 );
}

void MoonPosition::Reinitialize()
{
    this->longitudeIntegrated = 0.0;
    this->latitudeIntegrated = 0.0;
    this->earthToMoonRadius = 0.0;
    this->earthEccentricity = 0.0;
    this->actionOfVenus = 0.0;
    this->actionOfJupiter = 0.0;
    this->actionOfSomething = 0.0;
    this->moonMeanElongation = 0.0;
    this->moonMeanLongitude = 0.0;
    this->moonArgumentOfLatitude = 0.0;
    this->moonMeanAnomaly = 0.0;
    this->sunMeanAnomaly = 0.0;
    this->geocentricLongitude = 0.0;
    this->geocentricLatitude = 0.0;
    this->earthToMoonDistance = 0.0;
    this->cartesianCoordinates[ 0 ] = 0;
    this->cartesianCoordinates[ 1 ] = 0;
    this->cartesianCoordinates[ 2 ] = 0;
}

void MoonPosition::Iterate()
{
    CalcDynamicTime();
    CalcMoonMeanLongitude();
    CalcMoonMeanElongation();
    CalcSunMeanAnomaly();
    CalcMoonMeanAnomaly();
    CalcMoonArgumentOfLatitude();
    CalcActionOfVenus();
    CalcActionOfJupiter();
    CalcActionOfSomething();
    CalcEarthEccentricity();
}

void MoonPosition::Integrate()
{
    Reinitialize();
    Iterate();

    for( int i = 0; i < 60; ++i )
    {
        this->longitudeIntegrated += LONGITUDE_SINE_COEFF[i] *
                                     pow( earthEccentricity, floor( abs( LONGITUDE_PERIODIC_TERM_MATRIX[i][1] ) ) ) *
                                          sin( ( M_PI/180 ) * (
                                          LONGITUDE_PERIODIC_TERM_MATRIX[i][0] * moonMeanElongation +
                                          LONGITUDE_PERIODIC_TERM_MATRIX[i][1] * sunMeanAnomaly +
                                          LONGITUDE_PERIODIC_TERM_MATRIX[i][2] * moonMeanAnomaly +
                                          LONGITUDE_PERIODIC_TERM_MATRIX[i][3] * moonArgumentOfLatitude ) );

        this->earthToMoonRadius += LONGITUDE_COSINE_COEFF[i] *
                               pow( earthEccentricity, floor( abs( LONGITUDE_PERIODIC_TERM_MATRIX[i][1] ) ) ) *
                               cos( ( M_PI/180 ) * (
                                    LONGITUDE_PERIODIC_TERM_MATRIX[i][0] * moonMeanElongation +
                                    LONGITUDE_PERIODIC_TERM_MATRIX[i][1] * sunMeanAnomaly +
                                    LONGITUDE_PERIODIC_TERM_MATRIX[i][2] * moonMeanAnomaly +
                                    LONGITUDE_PERIODIC_TERM_MATRIX[i][3] * moonArgumentOfLatitude ) );

        this->latitudeIntegrated += LATITUDE_SINE_COEFF[i] *
                                    pow( earthEccentricity, floor( abs( LATITUDE_PERIODIC_TERM_MATRIX[i][1] ) ) ) *
                                    sin( ( M_PI/180 ) * (
                                        LATITUDE_PERIODIC_TERM_MATRIX[i][0] * moonMeanElongation +
                                        LATITUDE_PERIODIC_TERM_MATRIX[i][1] * sunMeanAnomaly +
                                        LATITUDE_PERIODIC_TERM_MATRIX[i][2] * moonMeanAnomaly +
                                        LATITUDE_PERIODIC_TERM_MATRIX[i][3] * moonArgumentOfLatitude ) );

    }

    this->longitudeIntegrated = longitudeIntegrated +
                                LONGITUDE_ADDITIVE_COEFFS[0] * sin( ( M_PI/180 ) * actionOfVenus ) +
                                LONGITUDE_ADDITIVE_COEFFS[1] * sin( ( M_PI/180 ) * ( moonMeanLongitude - moonArgumentOfLatitude ) ) +
                                LONGITUDE_ADDITIVE_COEFFS[2] * sin( ( M_PI/180 ) * actionOfJupiter );

    this->latitudeIntegrated = latitudeIntegrated -
                               LATITUDE_ADDITIVE_COEFFS[0] * sin( ( M_PI/180 ) * moonMeanLongitude ) + 
                               LATITUDE_ADDITIVE_COEFFS[1] * sin( ( M_PI/180 ) * actionOfSomething ) +
                               LATITUDE_ADDITIVE_COEFFS[2] * sin( ( M_PI/180 ) * ( actionOfVenus - moonArgumentOfLatitude ) ) +
                               LATITUDE_ADDITIVE_COEFFS[3] * sin( ( M_PI/180 ) * ( actionOfVenus + moonArgumentOfLatitude ) ) +
                               LATITUDE_ADDITIVE_COEFFS[4] * sin( ( M_PI/180 ) * ( moonMeanLongitude - moonMeanAnomaly ) ) -
                               LATITUDE_ADDITIVE_COEFFS[5] * sin( ( M_PI/180 ) * ( moonMeanLongitude + moonMeanAnomaly ) );

    this->geocentricLongitude = moonMeanLongitude + ( longitudeIntegrated / 1000000 );
    this->geocentricLatitude = ( latitudeIntegrated / 1000000 );
    this->earthToMoonDistance = 385000.56 + ( earthToMoonRadius / 1000 );

    this->cartesianCoordinates[0] = earthToMoonDistance *
                                    cos( ( M_PI/180 ) * geocentricLatitude ) *
                                    cos( ( M_PI/180 ) * geocentricLongitude );

    this->cartesianCoordinates[1] = earthToMoonDistance *
                                    cos( ( M_PI/180 ) * geocentricLatitude ) *
                                    sin( ( M_PI/180 ) * geocentricLongitude );

    this->cartesianCoordinates[2] = earthToMoonDistance * sin( ( M_PI/180 ) * geocentricLatitude );

    PrintCurrentPosition();
}

void MoonPosition::PrintCurrentPosition()
{
    int year;
    int month;
    double day;

    JD_to_Calendar_Date( this->julianEphemerisDay, &year, &month, &day );

    printf("Year: %u\n", year);
    printf("Month: %u\n", month);
    printf("Day: %f\n", day);
}

void MoonPosition::Shutdown()
{
    printf( "\n=============================================\n" );
    printf( "          Moon Position at Shutdown\n" );
    printf( "               X = %f\n",this->cartesianCoordinates[0] );
    printf( "               Y = %f\n",this->cartesianCoordinates[1] );
    printf( "               Z = %f\n",this->cartesianCoordinates[2] );
    printf( "=============================================\n\n" );
}
