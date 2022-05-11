#pragma once

#include "Math/Math.h"

class Grid;

class GridManager final
{
public:
	GridManager();

	~GridManager();

	void SetGrids(uint16 gridNum);
	Grid*** GetGrids();

	uint16 GetGridNum() const;

	void RestGrids();

private:
	void deleteGrids();

private:
	Grid*** mGrids = nullptr;
	uint16 mGridNum = 96;
};