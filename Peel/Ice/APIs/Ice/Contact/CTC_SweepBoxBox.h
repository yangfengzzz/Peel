///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Include Guard
#ifndef CTCSWEEPBOXBOX_H
#define CTCSWEEPBOXBOX_H

CONTACT_API bool SweepBoxBox(
        const OBB& box0, const OBB& box1, const Point& dir, float length, Point& hit, Point& normal, float& t);

#endif  // CTCSWEEPBOXBOX_H
