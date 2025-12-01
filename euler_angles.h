#pragma once

#include <type_traits>

namespace math
{

	enum euler_angle_type {
		// parity
		EULER_ANGLES_ZYXs,
		EULER_ANGLES_YZXs,
		// repetition
		EULER_ANGLES_XYXs,
		EULER_ANGLES_XZXs,
		// rotating
		EULER_ANGLES_XYZr,
		EULER_ANGLES_XZYr,
		EULER_ANGLES_XYXr,
		EULER_ANGLES_XZXr,

		EULER_ANGLES_XZYs,
		EULER_ANGLES_ZXYs,
		EULER_ANGLES_YZYs,
		EULER_ANGLES_YXYs,

		EULER_ANGLES_YZXr,
		EULER_ANGLES_YXZr,
		EULER_ANGLES_YZYr,
		EULER_ANGLES_YXYr,

		EULER_ANGLES_YXZs,
		EULER_ANGLES_XYZs,
		EULER_ANGLES_ZXZs,
		EULER_ANGLES_ZYZs,

		EULER_ANGLES_ZXYr,
		EULER_ANGLES_ZYXr,
		EULER_ANGLES_ZXZr,
		EULER_ANGLES_ZYZr,
	};

	constexpr euler_angle_type EULER_ANGLES_XYZ = EULER_ANGLES_XYZs;
	constexpr euler_angle_type EULER_ANGLES_XZY = EULER_ANGLES_XZYs;
	constexpr euler_angle_type EULER_ANGLES_YZX = EULER_ANGLES_YZXs;
	constexpr euler_angle_type EULER_ANGLES_YXZ = EULER_ANGLES_YXZs;
	constexpr euler_angle_type EULER_ANGLES_ZXY = EULER_ANGLES_ZXYs;
	constexpr euler_angle_type EULER_ANGLES_ZYX = EULER_ANGLES_ZYXs;

	namespace euler_angles
	{
		constexpr int get_parity(euler_angle_type type) noexcept
		{
			return type & 0x1;
		}
		constexpr bool get_repetition(euler_angle_type type) noexcept
		{
			return type & 0x2;
		}
		constexpr bool get_rotating(euler_angle_type type) noexcept
		{
			return type & 0x4;
		}

		constexpr int get_next_axis(euler_angle_type type, int i_axis) noexcept
		{
			return get_parity(type) ? (i_axis + 2) % 3 : (i_axis + 1) % 3;
		}
		constexpr int get_i_axis(euler_angle_type type) noexcept
		{
			return (type >> 3) & 0x3;
		}
		constexpr int get_j_axis(euler_angle_type type) noexcept
		{
			return get_next_axis(type, get_i_axis(type));
		}
		constexpr int get_k_axis(euler_angle_type type) noexcept
		{
			return get_next_axis(type, get_j_axis(type));
		}

		constexpr int get_axis0(euler_angle_type type) noexcept
		{
			return get_i_axis(type);
		}
		constexpr int get_axis1(euler_angle_type type) noexcept
		{
			return get_j_axis(type);
		}
		constexpr int get_axis2(euler_angle_type type) noexcept
		{
			if(get_repetition(type))
				return get_i_axis(type);
			else
				return get_k_axis(type);
		}
	}


	template<euler_angle_type type, typename T>
	struct teuler_angles {
		using elem_type = T;
		static constexpr euler_angle_type angle_type = type;
		T rot0, rot1, rot2;
	};


	namespace euler_angles
	{
		template<euler_angle_type type, typename T>
		constexpr
			std::enable_if_t<!get_repetition(type), T>
		heading(const teuler_angles<type, T>& euler) noexcept
		{
			if (get_rotating(type))
				return euler.rot0;
			else
				return euler.rot2;
		}
		template<euler_angle_type type, typename T>
		constexpr
			std::enable_if_t<!get_repetition(type), T>
		attitude(const teuler_angles<type, T>& euler) noexcept
		{
			return euler.rot1;
		}
		template<euler_angle_type type, typename T>
		constexpr
			std::enable_if_t<!get_repetition(type), T>
		bank(const teuler_angles<type, T>& euler) noexcept
		{
			if (get_rotating(type))
				return euler.rot2;
			else
				return euler.rot0;
		}

		template<euler_angle_type type, typename T>
		constexpr
			std::enable_if_t<!get_repetition(type), T>
		yaw(const teuler_angles<type, T>& euler) noexcept
		{
			if (get_rotating(type))
				return euler.rot0;
			else
				return euler.rot2;
		}
		template<euler_angle_type type, typename T>
		constexpr
			std::enable_if_t<!get_repetition(type), T>
		pitch(const teuler_angles<type, T>& euler) noexcept
		{
			return euler.rot1;
		}
		template<euler_angle_type type, typename T>
		constexpr
			std::enable_if_t<!get_repetition(type), T>
		roll(const teuler_angles<type, T>& euler) noexcept
		{
			if (get_rotating(type))
				return euler.rot2;
			else
				return euler.rot0;
		}


		template<euler_angle_type type, typename T>
		constexpr
			std::enable_if_t<!get_repetition(type), T>
		azimuth(const teuler_angles<type, T>& euler) noexcept
		{
			if (get_rotating(type))
				return euler.rot0;
			else
				return euler.rot2;
		}
		template<euler_angle_type type, typename T>
		constexpr
			std::enable_if_t<!get_repetition(type), T>
		elevation(const teuler_angles<type, T>& euler) noexcept
		{
			return euler.rot1;
		}
		template<euler_angle_type type, typename T>
		constexpr
			std::enable_if_t<!get_repetition(type), T>
		tilt(const teuler_angles<type, T>& euler) noexcept
		{
			if (get_rotating(type))
				return euler.rot2;
			else
				return euler.rot0;
		}

		template<euler_angle_type type, typename T>
		constexpr
			std::enable_if_t<get_repetition(type), T>
		precession(const teuler_angles<type, T>& euler) noexcept
		{
			return euler.rot2;
		}
		template<euler_angle_type type, typename T>
		constexpr
			std::enable_if_t<get_repetition(type), T>
		nutation(const teuler_angles<type, T>& euler) noexcept
		{
			return euler.rot1;
		}
		template<euler_angle_type type, typename T>
		constexpr
			std::enable_if_t<get_repetition(type), T>
		spin(const teuler_angles<type, T>& euler) noexcept
		{
			return euler.rot0;
		}
	}
}

