

// RooBarb
#include "XmlConfig.h"
#include "TaskEngine.h"
using namespace jdb;

// STL
#include <iostream>
#include <exception>

#include "ToyMaker/ToyEffMaker.h"
#include "JPsiEfficiency.h"
#include "PhiEfficiency.h"


#define LOGURU_IMPLEMENTATION 1
#include "vendor/loguru.h"

int main( int argc, char* argv[] ) {
	

	Logger::setGlobalLogLevel( "none" );

	
	TaskFactory::registerTaskRunner<ToyEffMaker>( "ToyEffMaker" );
	TaskFactory::registerTaskRunner<JPsiEfficiency>( "JPsiEfficiency" );
	TaskFactory::registerTaskRunner<PhiEfficiency>( "PhiEfficiency" );
	

	TaskEngine engine( argc, argv, "ToyEffMaker" );


	return 0;
}
