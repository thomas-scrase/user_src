

#ifndef TOOLS_HPP
#define TOOLS_HPP




void get_max( double &value, double &value_max) {

	if (value>value_max)
	{
		value_max = value;
	}
}



void get_min( double &value, double &value_min) {

	if (value<value_min)
	{
		value_min = value;
	}
}


#endif // TOOLS_HPP
