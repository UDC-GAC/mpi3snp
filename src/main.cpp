#include "SearchMI.h"

int main(int argc, char* argv[])
{

	double stime, etime;

	/*get the startup time*/
	stime = Utils::getSysTime();

	Options* options = new Options();

	/*parse the arguments*/
	if (!options->parse(argc, argv))
	{
		options->printUsage();
		return 0;
	}

	Search* search = new SearchMI(options);

	search->execute();

	etime = Utils::getSysTime();
	Utils::log("Overall time: %.2f seconds\n", etime - stime);

	delete search;
	delete options;
	return 0;
}
