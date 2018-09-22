

// RooBarb
#include "XmlConfig.h"
#include "TaskEngine.h"
using namespace jdb;

// STL
#include <iostream>
#include <exception>

#include "ToyMaker/ToyEffMaker.h"
#include "JPsiEfficiency.h"
#include "Psi2SEfficiency.h"
#include "PhiEfficiency.h"

#include "CocktailLikeEfficiency.h"
#include "VirtualPhotonEfficiency.h"


#define LOGURU_IMPLEMENTATION 1
#include "vendor/loguru.h"

int main( int argc, char* argv[] ) {
	

	Logger::setGlobalLogLevel( "none" );

	
	TaskFactory::registerTaskRunner<ToyEffMaker>( "ToyEffMaker" );
	TaskFactory::registerTaskRunner<JPsiEfficiency>( "JPsiEfficiency" );
	TaskFactory::registerTaskRunner<Psi2SEfficiency>( "Psi2SEfficiency" );
	TaskFactory::registerTaskRunner<PhiEfficiency>( "PhiEfficiency" );

	TaskFactory::registerTaskRunner<CocktailLikeEfficiency>( "CocktailLikeEfficiency" );
	TaskFactory::registerTaskRunner<VirtualPhotonEfficiency>( "VirtualPhotonEfficiency" );
	

	TaskEngine engine( argc, argv, "ToyEffMaker" );


	return 0;
}
