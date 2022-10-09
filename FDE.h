
template<typename T>
T Cal_one_sided_FDE(T step_size, T input, std::function<T(T)> f){
	return (f(input + step_size) - f(input))/step_size;
};

template<typename T>
T Cal_Central_FDE(T step_size, T input, std::function<T(T)> f){
	return (f(input + step_size) - f(input - step_size))/(2 * step_size);
};