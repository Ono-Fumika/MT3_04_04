#include <Novice.h>
#include <cmath>
#include <assert.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <imgui.h>

const char kWindowTitle[] = "LD2B_01_オノ_フミカ_MT3_";

// 3次元ベクトルの構造体
struct Vector3 { float x, y, z; };

// 4ｘ4行列の構造体
struct Matrix4x4 {
	float m[4][4];
};

// 球
struct Sphere {
	Vector3 center; // 中心点
	float radius; // 半径
};

// 平面
struct Plane {
	Vector3 normal; //法線
	float distance; //原点から法線までの距離
};

// ボール
struct Ball {
	// 位置
	Vector3 position;
	// 速度
	Vector3 velocity;
	// 加速度
	Vector3 acceleration;
	// 質量
	float mass;
	// 半径
	float radius;
	// 色
	unsigned int color;
};

// 加算
Vector3 Add(const Vector3& v1, const Vector3& v2) {
	Vector3 result = {};

	result.x = v1.x + v2.x;
	result.y = v1.y + v2.y;
	result.z = v1.z + v2.z;

	return result;
};
// 減産
Vector3 Subtract(const Vector3& v1, const Vector3& v2) {
	Vector3 result = {};

	result.x = v1.x - v2.x;
	result.y = v1.y - v2.y;
	result.z = v1.z - v2.z;

	return result;
};
// 減産
Vector3 Subtract(const Vector3& v, float scalar) {
	Vector3 result = {};

	result.x = v.x - scalar;
	result.y = v.y - scalar;
	result.z = v.z - scalar;

	return result;
}
// 割り算
Vector3 Division(float scalar, const Vector3& v) {
	Vector3 result = {};

	result.x = v.x / scalar;
	result.y = v.y / scalar;
	result.z = v.z / scalar;

	return result;
}
// スカラー倍
Vector3 Multiply(float scalar, const Vector3& v) {
	Vector3 result = {};

	result.x = v.x * scalar;
	result.y = v.y * scalar;
	result.z = v.z * scalar;

	return result;
};
// 内席
float Dot(const Vector3& v1, const Vector3& v2) {
	float result = {};
	result = (v1.x * v2.x) + (v1.y * v2.y) + (v1.z * v2.z);

	return result;
};
// 長さ（ノルム）
float Length(const Vector3& v) {
	float result = {};

	result = sqrtf(Dot(v, v));

	return result;
};
// 正規化
Vector3 Normalize(const Vector3& v) {
	Vector3 result = {};

	result.x = (v.x / Length(v));
	result.y = (v.y / Length(v));
	result.z = (v.z / Length(v));

	return result;
};
// 座標変換
Vector3 Transform(const Vector3& vector, const Matrix4x4& matrix) {
	// ｗ = 1 がデカルト座標系であるので（ｘ、ｙ、z,1）のベクトルとしてmatrixとの積をとる
	Vector3 result = {};

	result.x =
		vector.x * matrix.m[0][0] +
		vector.y * matrix.m[1][0] +
		vector.z * matrix.m[2][0] +
		1.0f * matrix.m[3][0];

	result.y =
		vector.x * matrix.m[0][1] +
		vector.y * matrix.m[1][1] +
		vector.z * matrix.m[2][1] +
		1.0f * matrix.m[3][1];

	result.z =
		vector.x * matrix.m[0][2] +
		vector.y * matrix.m[1][2] +
		vector.z * matrix.m[2][2] +
		1.0f * matrix.m[3][2];

	float w =
		vector.x * matrix.m[0][3] +
		vector.y * matrix.m[1][3] +
		vector.z * matrix.m[2][3] +
		1.0f * matrix.m[3][3];

	// ベクトルに対して基本的な操作を行う行列でｗが0になることは無い
	assert(w != 0.0f);
	// ｗ＝1がデカルト座標系であるので、ｗ除算することで同次座標をデカルト座標に戻す
	result.x /= w;
	result.y /= w;
	result.z /= w;

	return result;
};
// X軸回転行列
Matrix4x4 MakeRotateXMatrix(float radian) {
	Matrix4x4 result = {};

	result.m[0][0] = 1;
	result.m[1][1] = std::cos(radian);
	result.m[1][2] = std::sin(radian);
	result.m[2][1] = -std::sin(radian);
	result.m[2][2] = std::cos(radian);
	result.m[3][3] = 1;

	return result;
};
// Y軸回転行列
Matrix4x4 MakeRotateYMatrix(float radian) {
	Matrix4x4 result = {};

	result.m[0][0] = std::cos(radian);
	result.m[0][2] = -std::sin(radian);
	result.m[1][1] = 1;
	result.m[2][0] = std::sin(radian);
	result.m[2][2] = std::cos(radian);
	result.m[3][3] = 1;

	return result;
};
// Z軸回転行列
Matrix4x4 MakeRotateZMatrix(float radian) {
	Matrix4x4 result = {};

	result.m[0][0] = std::cos(radian);
	result.m[0][1] = std::sin(radian);
	result.m[1][0] = -std::sin(radian);;
	result.m[1][1] = std::cos(radian);
	result.m[2][2] = 1;
	result.m[3][3] = 1;

	return result;
};
// 行列の積
Matrix4x4 Multiply(const Matrix4x4& m1, const Matrix4x4& m2) {
	Matrix4x4 result = {};

	for (int x = 0; x < 4; x++) {
		for (int y = 0; y < 4; y++) {
			for (int a = 0; a < 4; a++) {
				result.m[x][y] += m1.m[x][a] * m2.m[a][y];
			}
		}
	}

	return result;
};
// アフィン変換
Matrix4x4 MakeAffineMatrix(const Vector3& scale, const Vector3& rotate, const Vector3& translate) {
	Matrix4x4 result = {};

	Matrix4x4 rotateXMatrix = MakeRotateXMatrix(rotate.x);
	Matrix4x4 rotateYMatrix = MakeRotateYMatrix(rotate.y);
	Matrix4x4 rotateZMatrix = MakeRotateZMatrix(rotate.z);
	Matrix4x4 rotateXYZMatrix = Multiply(rotateXMatrix, Multiply(rotateYMatrix, rotateZMatrix));


	result.m[0][0] = scale.x * rotateXYZMatrix.m[0][0];
	result.m[0][1] = scale.x * rotateXYZMatrix.m[0][1];
	result.m[0][2] = scale.x * rotateXYZMatrix.m[0][2];

	result.m[1][0] = scale.y * rotateXYZMatrix.m[1][0];
	result.m[1][1] = scale.y * rotateXYZMatrix.m[1][1];
	result.m[1][2] = scale.y * rotateXYZMatrix.m[1][2];

	result.m[2][0] = scale.z * rotateXYZMatrix.m[2][0];
	result.m[2][1] = scale.z * rotateXYZMatrix.m[2][1];
	result.m[2][2] = scale.z * rotateXYZMatrix.m[2][2];

	result.m[3][0] = translate.x;
	result.m[3][1] = translate.y;
	result.m[3][2] = translate.z;
	result.m[3][3] = 1;

	return result;
};
// 逆行列
Matrix4x4 Inverse(const Matrix4x4& m) {
	Matrix4x4 result = {};
	float determinant = 0;

	// 行列式 |A| を求める
	determinant =
		(m.m[0][0] * m.m[1][1] * m.m[2][2] * m.m[3][3]) + (m.m[0][0] * m.m[1][2] * m.m[2][3] * m.m[3][1]) + (m.m[0][0] * m.m[1][3] * m.m[2][1] * m.m[3][2]) -
		(m.m[0][0] * m.m[1][3] * m.m[2][2] * m.m[3][1]) - (m.m[0][0] * m.m[1][2] * m.m[2][1] * m.m[3][3]) - (m.m[0][0] * m.m[1][1] * m.m[2][3] * m.m[3][2]) -
		(m.m[0][1] * m.m[1][0] * m.m[2][2] * m.m[3][3]) - (m.m[0][2] * m.m[1][0] * m.m[2][3] * m.m[3][1]) - (m.m[0][3] * m.m[1][0] * m.m[2][1] * m.m[3][2]) +
		(m.m[0][3] * m.m[1][0] * m.m[2][2] * m.m[3][1]) + (m.m[0][2] * m.m[1][0] * m.m[2][1] * m.m[3][3]) + (m.m[0][1] * m.m[1][0] * m.m[2][3] * m.m[3][2]) +
		(m.m[0][1] * m.m[1][2] * m.m[2][0] * m.m[3][3]) + (m.m[0][2] * m.m[1][3] * m.m[2][0] * m.m[3][1]) + (m.m[0][3] * m.m[1][1] * m.m[2][0] * m.m[3][2]) -
		(m.m[0][3] * m.m[1][2] * m.m[2][0] * m.m[3][1]) - (m.m[0][2] * m.m[1][1] * m.m[2][0] * m.m[3][3]) - (m.m[0][1] * m.m[1][3] * m.m[2][0] * m.m[3][2]) -
		(m.m[0][1] * m.m[1][2] * m.m[2][3] * m.m[3][0]) - (m.m[0][2] * m.m[1][3] * m.m[2][1] * m.m[3][0]) - (m.m[0][3] * m.m[1][1] * m.m[2][2] * m.m[3][0]) +
		(m.m[0][3] * m.m[1][2] * m.m[2][1] * m.m[3][0]) + (m.m[0][2] * m.m[1][1] * m.m[2][3] * m.m[3][0]) + (m.m[0][1] * m.m[1][3] * m.m[2][2] * m.m[3][0]);

	// 逆行列を求める
	result.m[0][0] = (1 / determinant) *
		((m.m[1][1] * m.m[2][2] * m.m[3][3]) + (m.m[1][2] * m.m[2][3] * m.m[3][1]) + (m.m[1][3] * m.m[2][1] * m.m[3][2]) -
			(m.m[1][3] * m.m[2][2] * m.m[3][1]) - (m.m[1][2] * m.m[2][1] * m.m[3][3]) - (m.m[1][1] * m.m[2][3] * m.m[3][2]));
	result.m[0][1] = (1 / determinant) *
		(-(m.m[0][1] * m.m[2][2] * m.m[3][3]) - (m.m[0][2] * m.m[2][3] * m.m[3][1]) - (m.m[0][3] * m.m[2][1] * m.m[3][2]) +
			(m.m[0][3] * m.m[2][2] * m.m[3][1]) + (m.m[0][2] * m.m[2][1] * m.m[3][3]) + (m.m[0][1] * m.m[2][3] * m.m[3][2]));
	result.m[0][2] = (1 / determinant) *
		((m.m[0][1] * m.m[1][2] * m.m[3][3]) + (m.m[0][2] * m.m[1][3] * m.m[3][1]) + (m.m[0][3] * m.m[1][1] * m.m[3][2]) -
			(m.m[0][3] * m.m[1][2] * m.m[3][1]) - (m.m[0][2] * m.m[1][1] * m.m[3][3]) - (m.m[0][1] * m.m[1][3] * m.m[3][2]));
	result.m[0][3] = (1 / determinant) *
		(-(m.m[0][1] * m.m[1][2] * m.m[2][3]) - (m.m[0][2] * m.m[1][3] * m.m[2][1]) - (m.m[0][3] * m.m[1][1] * m.m[2][2]) +
			(m.m[0][3] * m.m[1][2] * m.m[2][1]) + (m.m[0][2] * m.m[1][1] * m.m[2][3]) + (m.m[0][1] * m.m[1][3] * m.m[2][2]));

	result.m[1][0] = (1 / determinant) *
		(-(m.m[1][0] * m.m[2][2] * m.m[3][3]) - (m.m[1][2] * m.m[2][3] * m.m[3][0]) - (m.m[1][3] * m.m[2][0] * m.m[3][2]) +
			(m.m[1][3] * m.m[2][2] * m.m[3][0]) + (m.m[1][2] * m.m[2][0] * m.m[3][3]) + (m.m[1][0] * m.m[2][3] * m.m[3][2]));
	result.m[1][1] = (1 / determinant) *
		((m.m[0][0] * m.m[2][2] * m.m[3][3]) + (m.m[0][2] * m.m[2][3] * m.m[3][0]) + (m.m[0][3] * m.m[2][0] * m.m[3][2]) -
			(m.m[0][3] * m.m[2][2] * m.m[3][0]) - (m.m[0][2] * m.m[2][0] * m.m[3][3]) - (m.m[0][0] * m.m[2][3] * m.m[3][2]));
	result.m[1][2] = (1 / determinant) *
		(-(m.m[0][0] * m.m[1][2] * m.m[3][3]) - (m.m[0][2] * m.m[1][3] * m.m[3][0]) - (m.m[0][3] * m.m[1][0] * m.m[3][2]) +
			(m.m[0][3] * m.m[1][2] * m.m[3][0]) + (m.m[0][2] * m.m[1][0] * m.m[3][3]) + (m.m[0][0] * m.m[1][3] * m.m[3][2]));
	result.m[1][3] = (1 / determinant) *
		((m.m[0][0] * m.m[1][2] * m.m[2][3]) + (m.m[0][2] * m.m[1][3] * m.m[2][0]) + (m.m[0][3] * m.m[1][0] * m.m[2][2]) -
			(m.m[0][3] * m.m[1][2] * m.m[2][0]) - (m.m[0][2] * m.m[1][0] * m.m[2][3]) - (m.m[0][0] * m.m[1][3] * m.m[2][2]));

	result.m[2][0] = (1 / determinant) *
		((m.m[1][0] * m.m[2][1] * m.m[3][3]) + (m.m[1][1] * m.m[2][3] * m.m[3][0]) + (m.m[1][3] * m.m[2][0] * m.m[3][1]) -
			(m.m[1][3] * m.m[2][1] * m.m[3][0]) - (m.m[1][1] * m.m[2][0] * m.m[3][3]) - (m.m[1][0] * m.m[2][3] * m.m[3][1]));
	result.m[2][1] = (1 / determinant) *
		(-(m.m[0][0] * m.m[2][1] * m.m[3][3]) - (m.m[0][1] * m.m[2][3] * m.m[3][0]) - (m.m[0][3] * m.m[2][0] * m.m[3][1]) +
			(m.m[0][3] * m.m[2][1] * m.m[3][0]) + (m.m[0][1] * m.m[2][0] * m.m[3][3]) + (m.m[0][0] * m.m[2][3] * m.m[3][1]));
	result.m[2][2] = (1 / determinant) *
		((m.m[0][0] * m.m[1][1] * m.m[3][3]) + (m.m[0][1] * m.m[1][3] * m.m[3][0]) + (m.m[0][3] * m.m[1][0] * m.m[3][1]) -
			(m.m[0][3] * m.m[1][1] * m.m[3][0]) - (m.m[0][1] * m.m[1][0] * m.m[3][3]) - (m.m[0][0] * m.m[1][3] * m.m[3][1]));
	result.m[2][3] = (1 / determinant) *
		(-(m.m[0][0] * m.m[1][1] * m.m[2][3]) - (m.m[0][1] * m.m[1][3] * m.m[2][0]) - (m.m[0][3] * m.m[1][0] * m.m[2][1]) +
			(m.m[0][3] * m.m[1][1] * m.m[2][0]) + (m.m[0][1] * m.m[1][0] * m.m[2][3]) + (m.m[0][0] * m.m[1][3] * m.m[2][1]));

	result.m[3][0] = (1 / determinant) *
		(-(m.m[1][0] * m.m[2][1] * m.m[3][2]) - (m.m[1][1] * m.m[2][2] * m.m[3][0]) - (m.m[1][2] * m.m[2][0] * m.m[3][1]) +
			(m.m[1][2] * m.m[2][1] * m.m[3][0]) + (m.m[1][1] * m.m[2][0] * m.m[3][2]) + (m.m[1][0] * m.m[2][2] * m.m[3][1]));
	result.m[3][1] = (1 / determinant) *
		((m.m[0][0] * m.m[2][1] * m.m[3][2]) + (m.m[0][1] * m.m[2][2] * m.m[3][0]) + (m.m[0][2] * m.m[2][0] * m.m[3][1]) -
			(m.m[0][2] * m.m[2][1] * m.m[3][0]) - (m.m[0][1] * m.m[2][0] * m.m[3][2]) - (m.m[0][0] * m.m[2][2] * m.m[3][1]));
	result.m[3][2] = (1 / determinant) *
		(-(m.m[0][0] * m.m[1][1] * m.m[3][2]) - (m.m[0][1] * m.m[1][2] * m.m[3][0]) - (m.m[0][2] * m.m[1][0] * m.m[3][1]) +
			(m.m[0][2] * m.m[1][1] * m.m[3][0]) + (m.m[0][1] * m.m[1][0] * m.m[3][2]) + (m.m[0][0] * m.m[1][2] * m.m[3][1]));
	result.m[3][3] = (1 / determinant) *
		((m.m[0][0] * m.m[1][1] * m.m[2][2]) + (m.m[0][1] * m.m[1][2] * m.m[2][0]) + (m.m[0][2] * m.m[1][0] * m.m[2][1]) -
			(m.m[0][2] * m.m[1][1] * m.m[2][0]) - (m.m[0][1] * m.m[1][0] * m.m[2][2]) - (m.m[0][0] * m.m[1][2] * m.m[2][1]));

	return result;
};
// 透視投影行列
Matrix4x4 MakePerspectiveFovMatrix(float fovY, float aspectRatio, float nearClip, float farClip) {
	Matrix4x4 result = {};
	float cot = 1 / (std::tan(fovY / 2));

	result.m[0][0] = (1 / aspectRatio) * cot;
	result.m[1][1] = cot;
	result.m[2][2] = farClip / (farClip - nearClip);
	result.m[2][3] = 1;
	result.m[3][2] = -nearClip * (farClip / (farClip - nearClip));

	return result;
};
// ビューポート変換行列
Matrix4x4 MakeViewportMatrix(float left, float top, float width, float height, float minDepth, float maxDepth) {
	Matrix4x4 result = {};

	result.m[0][0] = width / 2;
	result.m[1][1] = -(height / 2);
	result.m[2][2] = maxDepth - minDepth;
	result.m[3][0] = left + (width / 2);
	result.m[3][1] = top + (height / 2);
	result.m[3][2] = minDepth;
	result.m[3][3] = 1;

	return result;
};
// クロス積
Vector3 Cross(const Vector3& v1, const Vector3& v2) {
	Vector3 result = {};
	result.x = (v1.y * v2.z) - (v1.z * v2.y);
	result.y = (v1.z * v2.x) - (v1.x * v2.z);
	result.z = (v1.x * v2.y) - (v1.y * v2.x);

	return result;
}
// 正射影ベクトル
Vector3 Project(const Vector3& v1, const Vector3& v2) {
	Vector3 result;

	result = Multiply(Dot(v1, Normalize(v2)), Normalize(v2)); // (a*b^)b^

	return result;
}

// グリッド
void DrawGrid(const Matrix4x4& viewProjectionMatrix, const Matrix4x4& viewportMatrix) {
	const float kGridHalfWidth = 2.0f;// 半分
	const uint32_t kSubdivision = 10;// 分割数
	const float kGridEvery = (kGridHalfWidth * 2.0f) / float(kSubdivision);// 1つ分の長さ

	Vector3 zLineStart;
	Vector3 zLineEnd;
	Vector3 xLineStart;
	Vector3 xLineEnd;

	// 奥から手前の線を順々に引いていく
	for (uint32_t xIndex = 0; xIndex <= kSubdivision; ++xIndex) {
		// 視点と頂点を求める
		zLineStart = Vector3(xIndex * kGridEvery - kGridHalfWidth, 0.0f, kGridHalfWidth);
		zLineEnd = Vector3(xIndex * kGridEvery - kGridHalfWidth, 0.0f, -kGridHalfWidth);

		// 変換
		Matrix4x4 worldMatrixStart = MakeAffineMatrix({ 1.0f,1.0f,1.0f }, { 0.0f,0.0f,0.0f }, zLineStart);
		Matrix4x4 worldMatrixEnd = MakeAffineMatrix({ 1.0f,1.0f,1.0f }, { 0.0f,0.0f,0.0f }, zLineEnd);

		Matrix4x4 wvpMatrixStart = Multiply(worldMatrixStart, viewProjectionMatrix);
		Matrix4x4 wvpMatrixEnd = Multiply(worldMatrixEnd, viewProjectionMatrix);

		Vector3 nbcVertexStart = Transform(Vector3{}, wvpMatrixStart);
		Vector3 screenStartPoint = Transform(nbcVertexStart, viewportMatrix);
		Vector3 nbcVertexEnd = Transform(Vector3{}, wvpMatrixEnd);
		Vector3 screenEndPoint = Transform(nbcVertexEnd, viewportMatrix);

		// 描画
		if (xIndex == 5) {
			Novice::DrawLine((int)screenStartPoint.x, (int)screenStartPoint.y, (int)screenEndPoint.x, (int)screenEndPoint.y, BLACK);
		}
		else {
			Novice::DrawLine((int)screenStartPoint.x, (int)screenStartPoint.y, (int)screenEndPoint.x, (int)screenEndPoint.y, 0xAAAAAAFF);
		}
	}

	// 左から右
	for (uint32_t xIndex = 0; xIndex <= kSubdivision; ++xIndex) {
		// 視点と頂点を求める
		xLineStart = Vector3(kGridHalfWidth, 0.0f, xIndex * kGridEvery - kGridHalfWidth);
		xLineEnd = Vector3(-kGridHalfWidth, 0.0f, xIndex * kGridEvery - kGridHalfWidth);

		// 変換
		Matrix4x4 worldMatrixStart = MakeAffineMatrix({ 1.0f,1.0f,1.0f }, { 0.0f,0.0f,0.0f }, xLineStart);
		Matrix4x4 worldMatrixEnd = MakeAffineMatrix({ 1.0f,1.0f,1.0f }, { 0.0f,0.0f,0.0f }, xLineEnd);

		Matrix4x4 wvpMatrixStart = Multiply(worldMatrixStart, viewProjectionMatrix);
		Matrix4x4 wvpMatrixEnd = Multiply(worldMatrixEnd, viewProjectionMatrix);

		Vector3 nbcVertexStart = Transform(Vector3{}, wvpMatrixStart);
		Vector3 screenStartPoint = Transform(nbcVertexStart, viewportMatrix);
		Vector3 nbcVertexEnd = Transform(Vector3{}, wvpMatrixEnd);
		Vector3 screenEndPoint = Transform(nbcVertexEnd, viewportMatrix);

		// 描画
		if (xIndex == 5) {
			Novice::DrawLine((int)screenStartPoint.x, (int)screenStartPoint.y, (int)screenEndPoint.x, (int)screenEndPoint.y, BLACK);
		}
		else {
			Novice::DrawLine((int)screenStartPoint.x, (int)screenStartPoint.y, (int)screenEndPoint.x, (int)screenEndPoint.y, 0xAAAAAAFF);
		}
	}

}
// 球
void DrawSphere(const Sphere& sphere, const Matrix4x4& viewProjectionMatrix, const Matrix4x4& viewPortMatrix, uint32_t color) {
	float pi = float(M_PI);
	const uint32_t kSudivision = 10; //分割数
	const float kLatEvery = pi / kSudivision; //緯度分割一つ分の角度
	const float kLonEvery = (pi * 2) / kSudivision; //緯度分割二つ分の角度
	// world座標系でのa,b,cを求める
	Vector3 a, b, c;

	// 緯度の方向に分割
	for (uint32_t latindex = 0; latindex < kSudivision; ++latindex) {
		sphere;
		float lat = -pi / 2.0f + kLatEvery * latindex; //現在の緯度

		// 経度の方向に分割
		for (uint32_t lonIndex = 0; lonIndex < kSudivision; ++lonIndex) {
			float lon = lonIndex * kLonEvery; //現在の経度
			// ワールド座標系での頂点を求める
			a = {
				sphere.center.x + sphere.radius * (cosf(lat) * cosf(lon)),
				sphere.center.y + sphere.radius * sinf(lat),
				sphere.center.z + sphere.radius * (cosf(lat) * sinf(lon))
			};
			b = {
				sphere.center.x + sphere.radius * (cosf(lat + kLatEvery) * cosf(lon)),
				sphere.center.y + sphere.radius * sinf(lat + kLatEvery),
				sphere.center.z + sphere.radius * (cosf(lat + kLatEvery) * sin(lon))
			};
			c = {
				sphere.center.x + sphere.radius * (cosf(lat) * cosf(lon + kLonEvery)),
				sphere.center.y + sphere.radius * sinf(lat),
				sphere.center.z + sphere.radius * (cos(lat) * sinf(lon + kLonEvery))
			};

			// 変換
			Matrix4x4 worldMatrixA = MakeAffineMatrix({ 1.0f,1.0f,1.0f }, { 0.0f,0.0f,0.0f }, a);
			Matrix4x4 worldMatrixB = MakeAffineMatrix({ 1.0f,1.0f,1.0f }, { 0.0f,0.0f,0.0f }, b);
			Matrix4x4 worldMatrixC = MakeAffineMatrix({ 1.0f,1.0f,1.0f }, { 0.0f,0.0f,0.0f }, c);

			Matrix4x4 wvpMatrixA = Multiply(worldMatrixA, viewProjectionMatrix);
			Matrix4x4 wvpMatrixB = Multiply(worldMatrixB, viewProjectionMatrix);
			Matrix4x4 wvpMatrixC = Multiply(worldMatrixC, viewProjectionMatrix);

			Vector3 nbcVertexA = Transform(Vector3{}, wvpMatrixA);
			Vector3 nbcVertexB = Transform(Vector3{}, wvpMatrixB);
			Vector3 nbcVertexC = Transform(Vector3{}, wvpMatrixC);

			Vector3 screenA = Transform(nbcVertexA, viewPortMatrix);
			Vector3 screenB = Transform(nbcVertexB, viewPortMatrix);
			Vector3 screenC = Transform(nbcVertexC, viewPortMatrix);

			// 描画
			Novice::DrawLine((int)screenA.x, (int)screenA.y, (int)screenB.x, (int)screenB.y, color);
			Novice::DrawLine((int)screenA.x, (int)screenA.y, (int)screenC.x, (int)screenC.y, color);


		}
	}
}
// 平面と球の衝突判定
bool IsColision(const Sphere& sphere, const Plane& plane) {
	float distance;
	distance = Dot(plane.normal, sphere.center) - plane.distance;

	if (distance < 0) {
		distance *= -1;
	}

	if (distance <= sphere.radius) {
		return true;
	}
	else { return false; }
};
// 垂直なベクトルを求める
Vector3 Perpendicular(const Vector3& vector) {
	if (vector.x != 0.0f || vector.y != 0.0f) {
		return{ -vector.y,vector.x,0.0f };
	}
	return { 0.0f,-vector.z,vector.y };
}
// 平面の描画
void DrawPlane(const Plane& plane, const Matrix4x4& viewProjectionMatrix, const Matrix4x4& viewportMatrix, uint32_t color) {
	Vector3 center = Multiply(plane.distance, plane.normal);
	Vector3 perpendiculars[4];
	// 法線と垂直なベクトルを求める
	perpendiculars[0] = Normalize(Perpendicular(plane.normal));
	// 垂直なベクトルの逆ベクトルを求める
	perpendiculars[1] = { -perpendiculars[0].x,-perpendiculars[0].y,-perpendiculars[0].z };
	//  法線と垂直なベクトル
	perpendiculars[2] = Cross(plane.normal, perpendiculars[0]);
	// ↑の逆ベクトルを求める
	perpendiculars[3] = { -perpendiculars[2].x,-perpendiculars[2].y,-perpendiculars[2].z };

	Vector3 points[4];
	for (int32_t index = 0; index < 4; ++index) {
		Vector3 extend = Multiply(2.0f, perpendiculars[index]);
		Vector3 point = Add(center, extend);
		points[index] = Transform(Transform(point, viewProjectionMatrix), viewportMatrix);
	}

	// 描画
	Novice::DrawLine((int)points[0].x, (int)points[0].y, (int)points[2].x, (int)points[2].y, color);
	Novice::DrawLine((int)points[2].x, (int)points[2].y, (int)points[1].x, (int)points[1].y, color);
	Novice::DrawLine((int)points[1].x, (int)points[1].y, (int)points[3].x, (int)points[3].y, color);
	Novice::DrawLine((int)points[3].x, (int)points[3].y, (int)points[0].x, (int)points[0].y, color);
}

// 反射ベクトルを求める関数
Vector3 Reflect(const Vector3& input, const Vector3& nomal) {
	// 反射ベクトル
	Vector3 r;
	r = Subtract(input, Multiply(2 * (Dot(input, nomal)), nomal));
	return r;
}

// Windowsアプリでのエントリーポイント(main関数)
int WINAPI WinMain(HINSTANCE, HINSTANCE, LPSTR, int) {

	// ライブラリの初期化
	Novice::Initialize(kWindowTitle, 1280, 720);

	// キー入力結果を受け取る箱
	char keys[256] = {0};
	char preKeys[256] = {0};

	Vector3 cameraRotate{ 0.26f,0.0f,0.0f };
	Vector3 cameraTranslate{ 0.0f,1.9f,-6.49f };
	int kWindowWidth = 1280;
	int kWindowHeight = 720;

	// 平面
	Plane plane;
	plane.normal = Normalize({ -0.2f,0.9f,-0.3f });
	plane.distance = 0.0f;

	// ボール
	Ball ball{};
	ball.position = { 0.8f,1.2f,0.3f };
	ball.acceleration = { 0.0f,-9.8f,0.0f };
	ball.mass = 2.0f;
	ball.radius = 0.05f;
	ball.color = WHITE;

	// 球
	Sphere sphere = {};

	// 反発係数
	float e = 0.8f;

	// デルタタイム
	float deltaTime = 1.0f / 60.0f;

	// ウィンドウの×ボタンが押されるまでループ
	while (Novice::ProcessMessage() == 0) {
		// フレームの開始
		Novice::BeginFrame();

		// キー入力を受け取る
		memcpy(preKeys, keys, 256);
		Novice::GetHitKeyStateAll(keys);

		///
		/// ↓更新処理ここから
		///
		
		ball.velocity = Add(ball.velocity, Multiply(deltaTime, ball.acceleration));
		ball.position = Add(ball.position, Multiply(deltaTime, ball.velocity));
		if (IsColision(Sphere{ ball.position,ball.radius }, plane)) {
			Vector3 reflected = Reflect(ball.velocity, plane.normal);
			Vector3 projectToNomal = Project(reflected, plane.normal);
			Vector3 movingDirection = Subtract(reflected, projectToNomal);
			ball.velocity = Add(Multiply(e, projectToNomal), movingDirection);
		}

		
		sphere = { ball.position,0.05f };

		// カメラ
		Matrix4x4 cameraWorldMatrix = MakeAffineMatrix({ 1.0f,1.0f,1.0f }, cameraRotate, cameraTranslate);
		Matrix4x4 viewMatrix = Inverse(cameraWorldMatrix);
		Matrix4x4 projectionMatrix = MakePerspectiveFovMatrix(0.45f, float(kWindowWidth) / float(kWindowHeight), 0.1f, 100.0f);
		Matrix4x4 viewProjectionMatrix = Multiply(viewMatrix, projectionMatrix);
		Matrix4x4 viewportMatrix = MakeViewportMatrix(0, 0, float(kWindowWidth), float(kWindowHeight), 0.0f, 1.0f);

		///
		/// ↑更新処理ここまで
		///

		///
		/// ↓描画処理ここから
		///

		// グリッド描画
		DrawGrid(viewProjectionMatrix, viewportMatrix);
		// 振り子のボール
		DrawSphere(sphere, viewProjectionMatrix, viewportMatrix, WHITE);
		// 平面の描画
		DrawPlane(plane, viewProjectionMatrix, viewportMatrix, WHITE);

		///
		/// ↑描画処理ここまで
		///

		// フレームの終了
		Novice::EndFrame();

		// ESCキーが押されたらループを抜ける
		if (preKeys[DIK_ESCAPE] == 0 && keys[DIK_ESCAPE] != 0) {
			break;
		}
	}

	// ライブラリの終了
	Novice::Finalize();
	return 0;
}
