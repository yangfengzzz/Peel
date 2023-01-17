
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Include Guard
#ifndef ICESTACKTRACK_H
#define ICESTACKTRACK_H

#ifndef _WIN64
class ICECORE_API StackTrack {
public:
    StackTrack();
    udword GetUsedStackSize() const;

private:
    void* mStackPtr;
};
#endif

#endif  // ICESTACKTRACK_H
