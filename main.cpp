#include "euler_math.h"

int main()
{
	constexpr euler_angle_type type = EULER_ANGLES_ZXZs;
	using etype = double;
	teuler_angles<type, etype> euler{ 4.1, 0.2, 0.3 };
	tquaternion<etype> q1 = math::to_quaternion(euler);
	etype c1 = cos(euler.rot1 * 0.5);
	etype s1 = sin(euler.rot1 * 0.5);
	etype c0 = cos(euler.rot0 * 0.5);
	etype s0 = sin(euler.rot0 * 0.5);
	etype c2 = cos(euler.rot2 * 0.5);
	etype s2 = sin(euler.rot2 * 0.5);
	teuler_angles<type, etype> euler1 = math::to_euler_keep_sign<type>(q1);
	return 0;
}
