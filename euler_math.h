#pragma once

#include "euler_angles.h"
#include <cmath>


namespace math
{

	// use your own
	template<typename T>
	struct tvector3
	{
		T x, y, z;
		constexpr T operator[](size_t index) const noexcept
		{
			return *(&x + index);
		}
		T& operator[](size_t index) noexcept
		{
			return *(&x + index);
		}
	};

	// use your own
	template<typename T>
	struct tquaternion
	{
		T w, x, y, z;
	private:
		tquaternion(T w, T x, T y, T z) : w(w), x(x), y(y), z(z) {}
	public:
		tquaternion() = default;
		tquaternion(const tquaternion<T>&) = default;
		tquaternion(tquaternion<T>&&) = default;
		tquaternion<T>& operator=(const tquaternion<T>&) = default;
		tquaternion<T>& operator=(tquaternion<T>&&) = default;
		static constexpr tquaternion<T> from_real_pure(T real, const tvector3<T>& pure) noexcept
		{
			tquaternion<T> q;
			q.w = real;
			q.x = pure.x;
			q.y = pure.y;
			q.z = pure.z;
			return q;
		}

		static constexpr tquaternion<T>	wxyz(T w, T x, T y, T z) noexcept
		{
			return 	tquaternion<T>{w, x, y, z};
		}

		static constexpr tquaternion<T>	xyzw(T x, T y, T z, T w) noexcept
		{
			return 	tquaternion<T>{w, x, y, z};
		}

		constexpr T real() const noexcept
		{
			return w;
		}
		constexpr tvector3<T> pure() const noexcept
		{
			return tvector3<T>{ x, y, z };
		}
	};

	template<typename T>
	tquaternion<T> operator*(const tquaternion<T>& a, const tquaternion<T>& b) noexcept
	{
		return tquaternion<T>::wxyz(
			a.w * b.w - a.x * b.x - a.y * b.y - a.z * b.z,
			a.w * b.x + a.x * b.w + a.y * b.z - a.z * b.y,
			a.w * b.y + a.y * b.w + a.z * b.x - a.x * b.z,
			a.w * b.z + a.z * b.w + a.x * b.y - a.y * b.x);
	}

	// use your own, this one stores in row-major order
	template<typename T>
	struct tmatrix3
	{
	private:
		T m[3][3];
		constexpr tmatrix3(
			T m00, T m01, T m02,
			T m10, T m11, T m12,
			T m20, T m21, T m22) noexcept :
			m{
				{m00, m01, m02},
				{m10, m11, m12},
				{m20, m21, m22}
			}
		{
		}
	public:
		tmatrix3() = default;
		tmatrix3(const tmatrix3<T>&) = default;
		tmatrix3(tmatrix3<T>&&) = default;
		tmatrix3<T>& operator=(const tmatrix3<T>&) = default;
		tmatrix3<T>& operator=(tmatrix3<T>&&) = default;

		constexpr T operator()(size_t row, size_t column) const noexcept
		{
			return m[row][column];
		}
		T& operator()(size_t row, size_t column) noexcept
		{
			return m[row][column];
		}
		static tmatrix3<T> from_row_elements(
			T m00, T m01, T m02,
			T m10, T m11, T m12,
			T m20, T m21, T m22) noexcept
		{
			return tmatrix3<T>(
				m00, m01, m02,
				m10, m11, m12,
				m20, m21, m22
			);
		}
		static tmatrix3<T> from_column_elements(
			T m00, T m01, T m02,
			T m10, T m11, T m12,
			T m20, T m21, T m22) noexcept
		{
			return tmatrix3<T>(
				m00, m10, m20,
				m01, m11, m21,
				m02, m12, m22
			);
		}
	};

	template<typename T>
	tmatrix3<T> operator*(const tmatrix3<T>& a, const tmatrix3<T>& b) noexcept
	{
		return tmatrix3<T>::from_row_elements(
			a(0, 0) * b(0, 0) + a(0, 1) * b(1, 0) + a(0, 2) * b(2, 0),
			a(0, 0) * b(0, 1) + a(0, 1) * b(1, 1) + a(0, 2) * b(2, 1),
			a(0, 0) * b(0, 2) + a(0, 1) * b(1, 2) + a(0, 2) * b(2, 2),

			a(1, 0) * b(0, 0) + a(1, 1) * b(1, 0) + a(1, 2) * b(2, 0),
			a(1, 0) * b(0, 1) + a(1, 1) * b(1, 1) + a(1, 2) * b(2, 1),
			a(1, 0) * b(0, 2) + a(1, 1) * b(1, 2) + a(1, 2) * b(2, 2),

			a(2, 0) * b(0, 0) + a(2, 1) * b(1, 0) + a(2, 2) * b(2, 0),
			a(2, 0) * b(0, 1) + a(2, 1) * b(1, 1) + a(2, 2) * b(2, 1),
			a(2, 0) * b(0, 2) + a(2, 1) * b(1, 2) + a(2, 2) * b(2, 2)
		);
	}


	using std::sin;
	using std::cos;
	using std::atan2;
	using std::sqrt;
	using std::signbit;
	using std::abs;

	namespace constants
	{
		template<typename T>
		constexpr T pi() noexcept
		{
			return static_cast<T>(3.14159265358979323846);
		}

		template<typename T>
		constexpr T two_pi() noexcept
		{
			return static_cast<T>(6.28318530717958647692);
		}
		template<typename T>
		constexpr T half_pi() noexcept
		{
			return static_cast<T>(1.57079632679489661923);
		}
	}

	template<typename T>
	tmatrix3<T> identity_mat3() noexcept
	{
		tmatrix3<T> m;
		m(0, 0) = 1; m(0, 1) = 0; m(0, 2) = 0;
		m(1, 0) = 0; m(1, 1) = 1; m(1, 2) = 0;
		m(2, 0) = 0; m(2, 1) = 0; m(2, 2) = 1;
		return m;
	}

	template<typename T>
	tmatrix3<T> rotationX_mat3(T angle) noexcept
	{
		tmatrix3<T> m;
		T c = cos(angle);
		T s = sin(angle);
		m(0, 0) = 1; m(0, 1) = 0; m(0, 2) = 0;
		m(1, 0) = 0; m(1, 1) = c; m(1, 2) = -s;
		m(2, 0) = 0; m(2, 1) = s; m(2, 2) = c;
		return m;
	}
	template<typename T>
	tmatrix3<T> rotationY_mat3(T angle) noexcept
	{
		tmatrix3<T> m;
		T c = cos(angle);
		T s = sin(angle);
		m(0, 0) = c; m(0, 1) = 0; m(0, 2) = s;
		m(1, 0) = 0; m(1, 1) = 1; m(1, 2) = 0;
		m(2, 0) = -s; m(2, 1) = 0; m(2, 2) = c;
		return m;
	}

	template<typename T>
	tmatrix3<T> rotationZ_mat3(T angle) noexcept
	{
		tmatrix3<T> m;
		T c = cos(angle);
		T s = sin(angle);
		m(0, 0) = c; m(0, 1) = -s; m(0, 2) = 0;
		m(1, 0) = s; m(1, 1) = c; m(1, 2) = 0;
		m(2, 0) = 0; m(2, 1) = 0; m(2, 2) = 1;
		return m;
	}

	template<typename T>
	tmatrix3<T> rotation_i_mat3(size_t axis, T angle)
	{
		switch (axis)
		{
		case 0:
			return rotationX_mat3<T>(angle);
		case 1:
			return rotationY_mat3<T>(angle);
		case 2:
			return rotationZ_mat3<T>(angle);
		default:
			// we don't do platform specifics here
			//__debugbreak(); // invalid axis
			return identity_mat3<T>();
		}
	}
	template<typename T>
	tquaternion<T> identity_quat() noexcept
	{
		return tquaternion<T>::from_real_pure(static_cast<T>(1), tvector3<T>{ 0, 0, 0 });
	}
	template<typename T>
	tquaternion<T> rotationX_quat(T angle) noexcept
	{
		T r = static_cast<T>(0.5) * angle;
		return tquaternion<T>::from_real_pure(cos(r), tvector3<T>{ sin(r), 0, 0 });
	}
	template<typename T>
	tquaternion<T> rotationY_quat(T angle) noexcept
	{
		T r = static_cast<T>(0.5) * angle;
		return tquaternion<T>::from_real_pure(cos(r), tvector3<T>{ 0, sin(r), 0 });
	}
	template<typename T>
	tquaternion<T> rotationZ_quat(T angle) noexcept
	{
		T r = static_cast<T>(0.5) * angle;
		return tquaternion<T>::from_real_pure(cos(r), tvector3<T>{ 0, 0, sin(r) });
	}
	template<typename T>
	tquaternion<T> rotation_i_quat(size_t axis, T angle)
	{
		switch (axis)
		{
		case 0:
			return rotationX_quat<T>(angle);
		case 1:
			return rotationY_quat<T>(angle);
		case 2:
			return rotationZ_quat<T>(angle);
		default:
			// we don't do platform specifics here
			//__debugbreak(); // invalid axis
			return identity_quat<T>();
		}
	}

	template<typename T>
	constexpr tmatrix3<T> to_mat3(const tquaternion<T>& q) noexcept
	{
		tmatrix3<T> m;
		T xx = q.x * q.x;
		T yy = q.y * q.y;
		T zz = q.z * q.z;
		T xy = q.x * q.y;
		T xz = q.x * q.z;
		T yz = q.y * q.z;
		T wx = q.w * q.x;
		T wy = q.w * q.y;
		T wz = q.w * q.z;
		T ww = q.w * q.w;
		m(0, 0) = ww + xx - yy - zz;
		m(0, 1) = 2 * (xy - wz);
		m(0, 2) = 2 * (xz + wy);
		m(1, 0) = 2 * (xy + wz);
		m(1, 1) = ww - xx + yy - zz;
		m(1, 2) = 2 * (yz - wx);
		m(2, 0) = 2 * (xz - wy);
		m(2, 1) = 2 * (yz + wx);
		m(2, 2) = ww - xx - yy + zz;
		return m;
	}


	template<euler_angle_type type, typename T>
	tmatrix3<T> to_mat3(const teuler_angles<type, T>& euler) noexcept
	{
		constexpr int i = euler_angles::get_i_axis(type);
		constexpr int j = euler_angles::get_j_axis(type);
		constexpr int k = euler_angles::get_k_axis(type);
		static_assert(i >= 0 && j >= 0 && k >= 0);
		static_assert(i < 3 && j < 3 && k < 3);
		static_assert(i != j && j != k && i != k);
		T r0, r1, r2;
		if (euler_angles::get_rotating(type))
		{
			r0 = euler.rot2;
			r1 = euler.rot1;
			r2 = euler.rot0;
		}
		else
		{
			r0 = euler.rot0;
			r1 = euler.rot1;
			r2 = euler.rot2;
		}
		if (euler_angles::get_parity(type))
		{
			r0 = -r0;
			r1 = -r1;
			r2 = -r2;
		}
		const T c0 = cos(r0);
		const T s0 = sin(r0);
		const T c1 = cos(r1);
		const T s1 = sin(r1);
		const T c2 = cos(r2);
		const T s2 = sin(r2);
		const T cc = c0 * c2, cs = c0 * s2, sc = s0 * c2, ss = s0 * s2;
		if (euler_angles::get_repetition(type))
		{
			tmatrix3<T> m;
			m(i, i) = c1;		m(i, j) = s1 * s0;			m(i, k) = s1 * c0;
			m(j, i) = s1 * s2;	m(j, j) = -c1 * ss + cc;	m(j, k) = -c1 * cs - sc;
			m(k, i) = -s1 * c2;	m(k, j) = c1 * sc + cs;		m(k, k) = c1 * cc - ss;
			return m;
		}
		else
		{
			tmatrix3<T> m;
			m(i, i) = c1 * c2;	m(i, j) = s1 * sc - cs;	m(i, k) = s1 * cc + ss;
			m(j, i) = c1 * s2;	m(j, j) = s1 * ss + cc;	m(j, k) = s1 * cs - sc;
			m(k, i) = -s1;		m(k, j) = c1 * s0;		m(k, k) = c1 * c0;
			return m;
		}
	}

	template<euler_angle_type type, typename T, bool use_other = false>
	teuler_angles<type, T> to_euler_any(const tmatrix3<T>& m) noexcept
	{
		constexpr int i = euler_angles::get_i_axis(type);
		constexpr int j = euler_angles::get_j_axis(type);
		constexpr int k = euler_angles::get_k_axis(type);
		static_assert(i >= 0 && j >= 0 && k >= 0);
		static_assert(i < 3 && j < 3 && k < 3);
		static_assert(i != j && j != k && i != k);
		constexpr T p = euler_angles::get_parity(type) ? static_cast<T>(-1) : static_cast<T>(1);
		constexpr T sign_other = use_other ? static_cast<T>(-1) : static_cast<T>(1);
		T r0, r1, r2;
		if (euler_angles::get_repetition(type))
		{
			T s1 = sign_other * sqrt(m(j, i) * m(j, i) + m(k, i) * m(k, i));
			r1 = atan2(s1, m(i, i));
			r0 = p * atan2(sign_other * m(i, j), sign_other * m(i, k));
			T s0 = sin(r0); T c0 = cos(r0);
			r2 = p * atan2(sign_other * (c0 * m(k, j) - s0 * m(k, k)),
				sign_other * (c0 * m(j, j) - s0 * m(j, k)));
		}
		// XYZ type, has sine as the singular one
		else
		{
			T c1 = sign_other * sqrt(m(i, i) * m(i, i) + m(j, i) * m(j, i));
			r1 = p * atan2(-m(k, i), c1);
			r0 = p * atan2(sign_other * m(k, j), sign_other * m(k, k));
			T s0 = sin(r0); T c0 = cos(r0);
			r2 = p * atan2(sign_other * (-c0 * m(i, j) + s0 * m(i, k)),
				sign_other * (c0 * m(j, j) - s0 * m(j, k)));
		}
		teuler_angles<type, T> res;
		if (euler_angles::get_rotating(type))
		{
			res.rot0 = r2;
			res.rot1 = r1;
			res.rot2 = r0;
		}
		else
		{
			res.rot0 = r0;
			res.rot1 = r1;
			res.rot2 = r2;
		}
		return res;
	}

	template<euler_angle_type type, typename T>
	teuler_angles<type, T> to_euler(const tmatrix3<T>& m) noexcept
	{
		return to_euler_any<type, T, false>(m);
	}

	template<euler_angle_type type, typename T>
	teuler_angles<type, T> to_euler_other(const tmatrix3<T>& m) noexcept
	{
		return to_euler_any<type, T, true>(m);
	}

	template<euler_angle_type type, typename T>
	tquaternion<T> to_quaternion(const teuler_angles<type, T>& euler) noexcept
	{
		constexpr int i = euler_angles::get_i_axis(type);
		constexpr int j = euler_angles::get_j_axis(type);
		constexpr int k = euler_angles::get_k_axis(type);
		static_assert(i >= 0 && j >= 0 && k >= 0);
		static_assert(i < 3 && j < 3 && k < 3);
		static_assert(i != j && j != k && i != k);
		T r0, r1, r2;
		if (euler_angles::get_rotating(type))
		{
			r0 = static_cast<T>(0.5) * euler.rot2;
			r1 = static_cast<T>(0.5) * euler.rot1;
			r2 = static_cast<T>(0.5) * euler.rot0;
		}
		else
		{
			r0 = static_cast<T>(0.5) * euler.rot0;
			r1 = static_cast<T>(0.5) * euler.rot1;
			r2 = static_cast<T>(0.5) * euler.rot2;
		}
		if (euler_angles::get_parity(type))
		{
			r1 = -r1;
		}
		const T c1 = cos(r1);
		const T s1 = sin(r1);
		tvector3<T> pure;
		T real;
		if (euler_angles::get_repetition(type))
		{
			real = c1 * cos(r0 + r2);
			pure[i] = c1 * sin(r0 + r2);
			pure[j] = s1 * cos(r0 - r2);
			pure[k] = s1 * -sin(r0 - r2);
		}
		else
		{
			const T c0 = cos(r0);
			const T s0 = sin(r0);
			const T c2 = cos(r2);
			const T s2 = sin(r2);
			const T cc = c0 * c2, cs = c0 * s2, sc = s0 * c2, ss = s0 * s2;
			real = c1 * cc + s1 * ss;
			pure[i] = c1 * sc - s1 * cs;
			pure[j] = s1 * cc + c1 * ss;
			pure[k] = c1 * cs - s1 * sc;
		}
		if (euler_angles::get_parity(type))
		{
			pure[j] = -pure[j];
		}
		return tquaternion<T>::from_real_pure(real, pure);
	}

	template<typename T>
	T sign_from_bit(T x) noexcept
	{
		return signbit(x) ? static_cast<T>(-1) : static_cast<T>(1);
	}

	template<euler_angle_type type, typename T, bool use_other = false>
	teuler_angles<type, T> to_euler_keep_sign_any(const tquaternion<T>& q) noexcept
	{
		constexpr int i = euler_angles::get_i_axis(type);
		constexpr int j = euler_angles::get_j_axis(type);
		constexpr int k = euler_angles::get_k_axis(type);
		static_assert(i >= 0 && j >= 0 && k >= 0);
		static_assert(i < 3 && j < 3 && k < 3);
		static_assert(i != j && j != k && i != k);
		constexpr T p = euler_angles::get_parity(type) ? static_cast<T>(-1) : static_cast<T>(1);
		constexpr T sign_other = use_other ? static_cast<T>(-1) : static_cast<T>(1);
		T q0 = q.real();
		tvector3 <T> qv = q.pure();
		T q1 = qv[i];
		T q2 = qv[j];
		T q3 = qv[k];
		T r0, r1, r2;
		if (euler_angles::get_repetition(type))
		{
			T c1 = q0 * q0 + q1 * q1 - q2 * q2 - q3 * q3;
			T s_half1 = sqrt(q2 * q2 + q3 * q3);
			T c_half1 = sqrt(q0 * q0 + q1 * q1);
			r1 = sign_other * 2 * atan2(s_half1, c_half1);

			T half0p2, half0m2;
			half0p2 = atan2(q1, q0);
			half0m2 = p * atan2(-sign_other * q3, sign_other * q2);
			r0 = half0p2 + half0m2;
			r2 = half0p2 - half0m2;
		}
		else
		{
			/*
				q1*q3 - p*q0*q2 = -p*c1*s1
				q1*q3 + p*q0*q2 = 2*c0*s0*c2*s2+p*c1*s1*(c0^2-s0^2)(c2^2-s2^2)
				q0 + q2 = (c1+s1)(c0*c2+p*s0*s2)
				q0 - q2 = (c1-s1)(c0*c2-p*s0*s2)
				q1 + p *q3 = (c1-s1)(s0*c2 + p*c0*s2)
				q1 - p *q3 = (c1+s1)(s0*c2 - p*c0*s2)
				sin((r0+p*r2)/2) = (q1 + p *q3)/(c1-s1)
				cos((r0+p*r2)/2) = (q0 - q2)/(c1-s1)
				(r0+p*r2)/2 = atan2(q1 + p *q3, q0 - q2)
				(r0-p*r2)/2 = atan2(q1 - p *q3, q0 + q2)
			*/
			T s1_ = 2 * (q0 * q2 - p * q1 * q3);
			T f00 = q0 + q2, f01 = q0 - q2;
			T f10 = q1 - p * q3, f11 = q1 + p * q3;
			T c1_ = sign_other * sqrt((f11 * f11 + f01 * f01) * (f10 * f10 + f00 * f00));
			r1 = atan2(s1_, c1_);
			if (!use_other)
			{
				T half0p2 = atan2(f11, f01);
				T half0m2 = atan2(f10, f00);
				r0 = half0p2 + half0m2;
				r2 = p * (half0p2 - half0m2);
			}
			else
			{
				// r1 < -pi/2, (c1+s1) <= 0 (c1-s1) > 0
				if (signbit(s1_))
				{
					T half0p2 = atan2(f11, f01);
					T half0m2 = atan2(-f10, -f00);
					r0 = half0p2 + half0m2;
					r2 = p * (half0p2 - half0m2);
				}
				// r1 > pi/2, (c1+s1) > 0 (c1-s1) <= 0
				else
				{
					T half0p2 = atan2(-f11, -f01);
					T half0m2 = atan2(f10, f00);
					r0 = half0p2 + half0m2;
					r2 = p * (half0p2 - half0m2);
				}
			}
		}
		teuler_angles<type, T> res;
		if (euler_angles::get_rotating(type))
		{
			res.rot0 = r2;
			res.rot1 = r1;
			res.rot2 = r0;
		}
		else
		{
			res.rot0 = r0;
			res.rot1 = r1;
			res.rot2 = r2;
		}
		return res;

	}

	template<euler_angle_type type, typename T>
	teuler_angles<type, T> to_euler_keep_sign(const tquaternion<T>& q) noexcept
	{
		return to_euler_keep_sign_any<type, T, false>(q);
	}

	template<euler_angle_type type, typename T>
	teuler_angles<type, T> to_euler_keep_sign_other(const tquaternion<T>& q) noexcept
	{
		return to_euler_keep_sign_any<type, T, true>(q);
	}

	template<euler_angle_type type, typename T>
	teuler_angles<type, T> to_euler(const tquaternion<T>& q) noexcept
	{
		teuler_angles<type, T> su2_angles = to_euler_keep_sign<type>(q);
		if (su2_angles.rot0 < -constants::pi<T>())
		{
			su2_angles.rot0 += constants::two_pi();
		}
		else if (su2_angles.rot0 > constants::pi<T>())
		{
			su2_angles.rot0 -= constants::two_pi();
		}
		if (su2_angles.rot2 < -constants::pi<T>())
		{
			su2_angles.rot2 += constants::two_pi();
		}
		else if (su2_angles.rot2 > constants::pi<T>())
		{
			su2_angles.rot2 -= constants::two_pi();
		}
		return su2_angles;
	}

	template<euler_angle_type type, typename T>
	teuler_angles<type, T> to_euler_other(const tquaternion<T>& q) noexcept
	{
		teuler_angles<type, T> su2_angles = to_euler_keep_sign_other<type>(q);
		if (su2_angles.rot0 < -constants::pi<T>())
		{
			su2_angles.rot0 += constants::two_pi();
		}
		else if (su2_angles.rot0 > constants::pi<T>())
		{
			su2_angles.rot0 -= constants::two_pi();
		}
		if (su2_angles.rot2 < -constants::pi<T>())
		{
			su2_angles.rot2 += constants::two_pi();
		}
		else if (su2_angles.rot2 > constants::pi<T>())
		{
			su2_angles.rot2 -= constants::two_pi();
		}
		return su2_angles;
	}
}
