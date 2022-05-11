#include <cassert>

#include "Grid/GridManager.h"
#include "Grid/Grid.h"

GridManager::GridManager()
	: mGrids(nullptr)
{

}

GridManager::~GridManager()
{
	deleteGrids();
}

void GridManager::SetGrids(uint16 gridNum)
{
	assert(mGrids == nullptr, "Grids already allocated");
	mGridNum = gridNum;
	mGrids = new Grid**[mGridNum];
	for (size_t x = 0; x < mGridNum; x++)
	{
		mGrids[x] = new Grid*[mGridNum];
		for (size_t y = 0; y < mGridNum; y++)
		{
			mGrids[x][y] = new Grid[mGridNum];
			for (size_t z = 0; z < mGridNum; z++)
			{
				mGrids[x][y][z].SetSize(mGridNum);
			}
		}
	}
	assert(mGrids != nullptr);
}

Grid*** GridManager::GetGrids()
{
	assert(mGrids != nullptr);
	return mGrids;
}

uint16 GridManager::GetGridNum() const
{
	return mGridNum;
}

void GridManager::RestGrids()
{
	for (size_t x = 0; x < mGridNum; x++)
	{
		for (size_t y = 0; y < mGridNum; y++)
		{
			for (size_t z = 0; z < mGridNum; z++)
			{
				mGrids[x][y][z].Reset();
			}
		}
	}
}

void GridManager::deleteGrids()
{
	if (mGrids == nullptr)
	{
		return;
	}

	for (size_t x = 0; x < mGridNum; x++)
	{
		for (size_t y = 0; y < mGridNum; y++)
		{
			if (mGrids[x][y] != nullptr)
			{
				delete[] mGrids[x][y];
			}
		}
		if (mGrids[x] != nullptr)
		{
			delete[] mGrids[x];
		}
	}
	delete[] mGrids;
	mGridNum = 0;
}