#include "euler_math.h"


using namespace math;
int main()
{
	constexpr euler_angle_type type = EULER_ANGLES_ZYZs;
	using etype = double;
	teuler_angles<type, etype> e0{ 0.1, -3.2, 0.3 };
	tquaternion<etype> q0 = to_quaternion(e0);
	teuler_angles<type, etype> e1 = to_euler_SU2<type, 3>(q0);
	tmatrix3<etype> m0 = to_mat3(e0);
	teuler_angles<type, etype> e2 = to_euler<type>(m0);
	return 0;
}
