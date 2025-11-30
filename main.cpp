#include "euler_math.h"


using namespace math;
int main()
{
	constexpr euler_angle_type type = EULER_ANGLES_ZXY;
	using etype = double;
	teuler_angles<type, etype> euler{ 4.1, 0.5, 0.3 };
	tquaternion<etype> q1 = to_quaternion(euler);
	teuler_angles<type, etype> euler1 = to_euler_keep_sign<type>(q1);
	return 0;
}
