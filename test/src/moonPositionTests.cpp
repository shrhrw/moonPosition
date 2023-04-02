#include <gtest/gtest.h>
#include "moonPosition.h"

class MoonPositionTest : public ::testing::Test {
    protected:
        MoonPosition moon;
};

/*
 * Initialization test, with all state variables set to 0.0 to test model set up.
 */
TEST_F( MoonPositionTest, InitializationTest ) {
    moon.Reinitialize();
    EXPECT_NEAR( 0.0, moon.longitudeIntegrated, 1e-6 );
    EXPECT_NEAR( 0.0, moon.latitudeIntegrated, 1e-6 );
    EXPECT_NEAR( 0.0, moon.earthToMoonRadius, 1e-6 );
    EXPECT_NEAR( 0.0, moon.dynamicalTime, 1e-6 );
    EXPECT_NEAR( 0.0, moon.earthEccentricity, 1e-6 );
    EXPECT_NEAR( 0.0, moon.actionOfVenus, 1e-6 );
    EXPECT_NEAR( 0.0, moon.actionOfJupiter, 1e-6 );
    EXPECT_NEAR( 0.0, moon.actionOfSomething, 1e-6 );
    EXPECT_NEAR( 0.0, moon.moonMeanElongation, 1e-6 );
    EXPECT_NEAR( 0.0, moon.moonMeanLongitude, 1e-6 );
    EXPECT_NEAR( 0.0, moon.moonArgumentOfLatitude, 1e-6 );
    EXPECT_NEAR( 0.0, moon.moonMeanAnomaly, 1e-6 );
    EXPECT_NEAR( 0.0, moon.sunMeanAnomaly, 1e-6 );
    EXPECT_NEAR( 0.0, moon.geocentricLongitude, 1e-6 );
    EXPECT_NEAR( 0.0, moon.geocentricLatitude, 1e-6 );
    EXPECT_NEAR( 0.0, moon.earthToMoonDistance, 1e-6 );
    EXPECT_NEAR( 0.0, moon.julianEphemerisDay, 1e-6 );
    EXPECT_NEAR( 0.0, moon.cartesianCoordinates[0], 1e-6 );
    EXPECT_NEAR( 0.0, moon.cartesianCoordinates[1], 1e-6 );
    EXPECT_NEAR( 0.0, moon.cartesianCoordinates[2], 1e-6 );
}

/*
 * Moon Mean Elongation function test.
 */
TEST_F( MoonPositionTest, MoonMeanElongationTest ) {
    moon.julianEphemerisDay = 2448724.5;
    moon.CalcDynamicTime();
    moon.CalcMoonMeanElongation();
    EXPECT_NEAR( 113.842, moon.moonMeanElongation, 0.001 );
}

/*
 * Sun Mean Anomaly function test.
 */
TEST_F( MoonPositionTest, SunMeanAnomalyTest ) {
    moon.julianEphemerisDay = 2448724.5;
    moon.CalcDynamicTime();
    moon.CalcSunMeanAnomaly();
    EXPECT_NEAR( 97.643514, moon.sunMeanAnomaly, 1e-6 );
}

/*
 * Moon Mean Anomaly function test.
 */
TEST_F( MoonPositionTest, MoonMeanAnomalyTest ) {
    moon.julianEphemerisDay = 2448724.5;
    moon.CalcDynamicTime();
    moon.CalcMoonMeanAnomaly();
    EXPECT_NEAR( 5.150833, moon.moonMeanAnomaly, 1e-6 );
}

/*
 * Moon Argument of Latitude function test.
 */
TEST_F( MoonPositionTest, MoonArgumentOfLatitudeTest ) {
    moon.julianEphemerisDay = 2448724.5;
    moon.CalcDynamicTime();
    moon.CalcMoonArgumentOfLatitude();
    EXPECT_NEAR( 219.889721, moon.moonArgumentOfLatitude, 1e-6 );
}

/*
 * Action of Venus function test.
 */
TEST_F( MoonPositionTest, ActionOfVenusTest ) {
    moon.julianEphemerisDay = 2448724.5;
    moon.CalcDynamicTime();
    moon.CalcActionOfVenus();
    EXPECT_NEAR( 109.57, moon.actionOfVenus, 0.01 );
}

/*
 * Action of Juptiter function test.
 */
TEST_F( MoonPositionTest, ActionOfJupiterTest ) {
    moon.julianEphemerisDay = 2448724.5;
    moon.CalcDynamicTime();
    moon.CalcActionOfJupiter();
    EXPECT_NEAR( 123.78, moon.actionOfJupiter, 0.01 );
}

/*
 * Action of Something function test.
 */
TEST_F( MoonPositionTest, ActionOfSomethingTest ) {
    moon.julianEphemerisDay = 2448724.5;
    moon.CalcDynamicTime();
    moon.CalcActionOfSomething();
    EXPECT_NEAR( 229.53, moon.actionOfSomething, 0.01 );
}

/*
 * Earth Eccentricity function test.
 */
TEST_F( MoonPositionTest, EarthEccentricityTest ) {
    moon.julianEphemerisDay = 2448724.5;
    moon.CalcDynamicTime();
    moon.CalcEarthEccentricity();
    EXPECT_NEAR( 1.000194, moon.earthEccentricity, 1e-6 );
}

/*
 * Full integration test, using Julian Ephemeris Day outlined in Meeus, page 342
 */
TEST_F( MoonPositionTest, IntegrateTest ) {
    moon.julianEphemerisDay = 2448724.5;
    moon.Iterate();
    moon.Integrate();
    EXPECT_NEAR( -1127527.033355, moon.longitudeIntegrated, 0.1 );
    EXPECT_NEAR( -3224463.902493, moon.latitudeIntegrated, 0.1 );
    EXPECT_NEAR( -16590875.183568, moon.earthToMoonRadius, 0.1 );
    EXPECT_NEAR( 133.162655, moon.geocentricLongitude, 1e-6 );
    EXPECT_NEAR( -3.224464, moon.geocentricLatitude, 1e-6 );
    EXPECT_NEAR( 368409.684816, moon.earthToMoonDistance, 1e-6 );
    EXPECT_NEAR( -251619.697297, moon.cartesianCoordinates[ 0 ], 1e-6 );
    EXPECT_NEAR( 268297.992270, moon.cartesianCoordinates[ 1 ], 1e-6 );
    EXPECT_NEAR( -20722.237874, moon.cartesianCoordinates[ 2 ], 1e-6 );
}

/*
 * Reduce functional tests (parameterized).
 */
class MoonPositionReduceTests : public ::testing::TestWithParam<std::tuple< double, double > > {
    protected:
        MoonPositionReduceTests() {}
        ~MoonPositionReduceTests() {}

        void SetUp() override {}

        void TearDown() override {}
            
        MoonPosition moon;
        double testValue = std::get<0>( GetParam() );
        double expectedValue = std::get<1>( GetParam() );
};

TEST_P( MoonPositionReduceTests, ReduceTests ) {
    EXPECT_NEAR( expectedValue, moon.Reduce( testValue ), 1e-6 );
}

INSTANTIATE_TEST_SUITE_P(
    ReduceTests,
    MoonPositionReduceTests,
    ::testing::Values(
        std::make_tuple( 721.3, 1.3 ),  // Value > 2pi (currently has rounding error)
        std::make_tuple( 13.6, 13.6),   // 0 < Value < 2pi 
        std::make_tuple( -13.6, 346.4), // Value < 0
        std::make_tuple( 0.0, 0.0),     // Value == 0
        std::make_tuple( 360.0, 0.0 )   // Value == 2pi
        ));

/*
 * Dynamic Time functional tests (parameterized).
 */
class MoonPositionDynamicTimeTests : public ::testing::TestWithParam<std::tuple< double, double > > {
    protected:
        MoonPositionDynamicTimeTests() {}
        ~MoonPositionDynamicTimeTests() {}

        void SetUp() override {}

        void TearDown() override {}
            
        MoonPosition moon;
        double testValue = std::get<0>( GetParam() );
        double expectedValue = std::get<1>( GetParam() );
};

TEST_P( MoonPositionDynamicTimeTests, DynamicTimeTests ) {
    moon.julianEphemerisDay = testValue;
    moon.CalcDynamicTime();
    EXPECT_NEAR( expectedValue, moon.dynamicalTime, 1e-6 );
}

INSTANTIATE_TEST_SUITE_P(
    DynamicTimeTests,
    MoonPositionDynamicTimeTests,
    ::testing::Values(
        std::make_tuple( 2448724.5, -0.0772211 )  
        ));

int main( int argc, char **argv ) {
    ::testing::InitGoogleTest( &argc, argv );
    return RUN_ALL_TESTS();
}
