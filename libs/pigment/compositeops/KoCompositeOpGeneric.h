/*
 * SPDX-FileCopyrightText: 2011 Silvio Heinrich <plassy@web.de>
 *
 * SPDX-License-Identifier: LGPL-2.1-or-later
*/
#ifndef _KOCOMPOSITEOP_GENERIC_H_
#define _KOCOMPOSITEOP_GENERIC_H_

#include "KoCompositeOpFunctions.h"
#include "KoCompositeOpBase.h"
#include "KoCompositeOpGenericFunctorBase.h"

namespace detail {

/**
 * A special class to convert old-style composite function into a functor
 * with clamping properties
 */

template <class Traits,
          typename Traits::channels_type compositeFunc(typename Traits::channels_type, typename Traits::channels_type)>
struct CompositeFunctionWrapper : KoCompositeOpGenericFunctorBase<typename Traits::channels_type>
{
    using channels_type = typename Traits::channels_type;

    static inline channels_type composeChannel(channels_type src, channels_type dst) {
        return compositeFunc(src, dst);
    }
};

template <class Traits,
          void compositeFunc(float, float, float, float&, float&, float&)>
struct CompositeFunctionWrapperHSL : KoCompositeOpGenericFunctorBase<typename Traits::channels_type>
{
    using channels_type = typename Traits::channels_type;

    static inline void composeChannels(float sr, float sg, float sb, float &dr, float &dg, float &db) {
        return compositeFunc(sr, sg, sb, dr, dg, db);
    }
};

}

/**
 * Generic CompositeOp for separable channel compositing functions
 *
 * A template to generate a KoCompositeOp class by just specifying a
 * blending/compositing function. This template works with compositing functions
 * for separable channels (means each channel of a pixel can be processed separately)
 */
template<
    class Traits,
    typename CompositeOpFunctor,
    typename BlendingPolicy
>
class KoCompositeOpGenericSCFunctor: public KoCompositeOpBase< Traits, KoCompositeOpGenericSCFunctor<Traits,CompositeOpFunctor,BlendingPolicy> >
{
    typedef KoCompositeOpBase< Traits, KoCompositeOpGenericSCFunctor<Traits,CompositeOpFunctor,BlendingPolicy> > base_class;
    typedef typename Traits::channels_type                                            channels_type;

    static const qint32 channels_nb = Traits::channels_nb;
    static const qint32 alpha_pos   = Traits::alpha_pos;

public:
    KoCompositeOpGenericSCFunctor(const KoColorSpace* cs, const QString& id, const QString& category)
        : base_class(cs, id, category) { }

public:
    template<bool alphaLocked, bool allChannelFlags>
    inline static channels_type composeColorChannels(const channels_type* src, channels_type srcAlpha,
                                                     channels_type*       dst, channels_type dstAlpha, channels_type maskAlpha,
                                                     channels_type opacity, const QBitArray& channelFlags) {
        using namespace Arithmetic;

        srcAlpha = mul(srcAlpha, maskAlpha, opacity);

        if (isZeroValueFuzzy(srcAlpha)) {
            return dstAlpha;
        }

        if(alphaLocked) {
            if(!isZeroValueFuzzy(dstAlpha)) {
                for(qint32 i=0; i <channels_nb; i++) {
                    if(i != alpha_pos && (allChannelFlags || channelFlags.testBit(i))) {
                        const channels_type srcInBlendSpace =
                                CompositeOpFunctor::clampSourceChannelValue(
                                    BlendingPolicy::toAdditiveSpace(
                                        src[i]));
                        const channels_type dstInBlendSpace =
                                CompositeOpFunctor::clampDestinationChannelValue(
                                    BlendingPolicy::toAdditiveSpace(
                                        dst[i]));

                        dst[i] = BlendingPolicy::fromAdditiveSpace(
                            lerp(dstInBlendSpace,
                                 CompositeOpFunctor::composeChannel(srcInBlendSpace, dstInBlendSpace),
                                 srcAlpha));
                    }
                }
            }

            return dstAlpha;
        } else if (isZeroValueFuzzy(dstAlpha)) {
            for(qint32 i=0; i <channels_nb; i++) {
                if(i != alpha_pos && (allChannelFlags || channelFlags.testBit(i))) {
                    dst[i] = BlendingPolicy::fromAdditiveSpace(
                                CompositeOpFunctor::clampSourceChannelValue(
                                    BlendingPolicy::toAdditiveSpace(src[i])));
                }
            }
            return srcAlpha;
        } else if (isUnitValueFuzzy(dstAlpha)) {
            for(qint32 i=0; i <channels_nb; i++) {
                if(i != alpha_pos && (allChannelFlags || channelFlags.testBit(i))) {
                    const channels_type srcInBlendSpace =
                            CompositeOpFunctor::clampSourceChannelValue(
                                BlendingPolicy::toAdditiveSpace(
                                    src[i]));
                    const channels_type dstInBlendSpace =
                            CompositeOpFunctor::clampDestinationChannelValue(
                                BlendingPolicy::toAdditiveSpace(
                                    dst[i]));

                    dst[i] = BlendingPolicy::fromAdditiveSpace(
                        lerp(dstInBlendSpace,
                             CompositeOpFunctor::composeChannel(srcInBlendSpace, dstInBlendSpace),
                             srcAlpha));
                }
            }
            return unitValue<channels_type>();
        }  else if (isUnitValueFuzzy(srcAlpha)) {
            for(qint32 i=0; i <channels_nb; i++) {
                if(i != alpha_pos && (allChannelFlags || channelFlags.testBit(i))) {
                    const channels_type srcInBlendSpace =
                            CompositeOpFunctor::clampSourceChannelValue(
                                BlendingPolicy::toAdditiveSpace(
                                    src[i]));
                    const channels_type dstInBlendSpace =
                            CompositeOpFunctor::clampDestinationChannelValue(
                                BlendingPolicy::toAdditiveSpace(
                                    dst[i]));

                    dst[i] = BlendingPolicy::fromAdditiveSpace(
                        lerp(srcInBlendSpace,
                             CompositeOpFunctor::composeChannel(srcInBlendSpace, dstInBlendSpace),
                             dstAlpha));
                }
            }
            return unitValue<channels_type>();
        } else {
            channels_type newDstAlpha = unionShapeOpacity(srcAlpha, dstAlpha);

            if (!isZeroValueFuzzy(newDstAlpha)) {

                for(qint32 i=0; i <channels_nb; i++) {
                    if(i != alpha_pos && (allChannelFlags || channelFlags.testBit(i))) {
                        const channels_type srcInBlendSpace =
                                CompositeOpFunctor::clampSourceChannelValue(
                                    BlendingPolicy::toAdditiveSpace(
                                        src[i]));
                        const channels_type dstInBlendSpace =
                                CompositeOpFunctor::clampDestinationChannelValue(
                                    BlendingPolicy::toAdditiveSpace(
                                        dst[i]));

                        channels_type result =
                            blend(srcInBlendSpace, srcAlpha,
                                  dstInBlendSpace, dstAlpha,
                                  CompositeOpFunctor::composeChannel(srcInBlendSpace, dstInBlendSpace));

                        dst[i] = BlendingPolicy::fromAdditiveSpace(div(result, newDstAlpha));
                    }
                }
            }

            return newDstAlpha;
        }
    }
};

template<
    class Traits,
    typename Traits::channels_type compositeFunc(typename Traits::channels_type, typename Traits::channels_type),
    typename BlendingPolicy
>
class KoCompositeOpGenericSC : public KoCompositeOpGenericSCFunctor<Traits, detail::CompositeFunctionWrapper<Traits, compositeFunc>, BlendingPolicy>
{
protected:
    using base_class = KoCompositeOpGenericSCFunctor<Traits, detail::CompositeFunctionWrapper<Traits, compositeFunc>, BlendingPolicy>;
public:
    using base_class::base_class;
};

/**
 * Generic CompositeOp for nonseparable/HSL channel compositing functions
 *
 * A template to generate a KoCompositeOp class by just specifying a
 * blending/compositing function. This template works with compositing functions
 * for RGB channels only (the channels can not be processed separately)
 */
template<class Traits, typename CompositeOpFunctor>
class KoCompositeOpGenericHSLFunctor: public KoCompositeOpBase< Traits, KoCompositeOpGenericHSLFunctor<Traits,CompositeOpFunctor> >
{
    typedef KoCompositeOpBase< Traits, KoCompositeOpGenericHSLFunctor<Traits,CompositeOpFunctor>> base_class;
    typedef typename Traits::channels_type channels_type;

    static const qint32 red_pos   = Traits::red_pos;
    static const qint32 green_pos = Traits::green_pos;
    static const qint32 blue_pos  = Traits::blue_pos;

public:
    KoCompositeOpGenericHSLFunctor(const KoColorSpace* cs, const QString& id, const QString& category)
        : base_class(cs, id, category) { }

public:
    template<bool alphaLocked, bool allChannelFlags>
    inline static channels_type composeColorChannels(const channels_type* src, channels_type srcAlpha,
                                                     channels_type*       dst, channels_type dstAlpha, channels_type maskAlpha,
                                                     channels_type opacity, const QBitArray& channelFlags) {
        using namespace Arithmetic;

        srcAlpha = mul(srcAlpha, maskAlpha, opacity);

        if(alphaLocked) {
            if(dstAlpha != zeroValue<channels_type>()) {
                float srcR = scale<float>(CompositeOpFunctor::clampSourceChannelValue(src[red_pos]));
                float srcG = scale<float>(CompositeOpFunctor::clampSourceChannelValue(src[green_pos]));
                float srcB = scale<float>(CompositeOpFunctor::clampSourceChannelValue(src[blue_pos]));

                float dstR = scale<float>(CompositeOpFunctor::clampDestinationChannelValue(dst[red_pos]));
                float dstG = scale<float>(CompositeOpFunctor::clampDestinationChannelValue(dst[green_pos]));
                float dstB = scale<float>(CompositeOpFunctor::clampDestinationChannelValue(dst[blue_pos]));

                CompositeOpFunctor::composeChannels(srcR, srcG, srcB, dstR, dstG, dstB);

                if(allChannelFlags || channelFlags.testBit(red_pos))
                    dst[red_pos] = lerp(dst[red_pos], scale<channels_type>(dstR), srcAlpha);

                if(allChannelFlags || channelFlags.testBit(green_pos))
                    dst[green_pos] = lerp(dst[green_pos], scale<channels_type>(dstG), srcAlpha);

                if(allChannelFlags || channelFlags.testBit(blue_pos))
                    dst[blue_pos] = lerp(dst[blue_pos], scale<channels_type>(dstB), srcAlpha);
            }

            return dstAlpha;
        }
        else {
            channels_type newDstAlpha = unionShapeOpacity(srcAlpha, dstAlpha);

            if(newDstAlpha != zeroValue<channels_type>()) {
                float srcR = scale<float>(CompositeOpFunctor::clampSourceChannelValue(src[red_pos]));
                float srcG = scale<float>(CompositeOpFunctor::clampSourceChannelValue(src[green_pos]));
                float srcB = scale<float>(CompositeOpFunctor::clampSourceChannelValue(src[blue_pos]));

                float dstR = scale<float>(CompositeOpFunctor::clampDestinationChannelValue(dst[red_pos]));
                float dstG = scale<float>(CompositeOpFunctor::clampDestinationChannelValue(dst[green_pos]));
                float dstB = scale<float>(CompositeOpFunctor::clampDestinationChannelValue(dst[blue_pos]));

                CompositeOpFunctor::composeChannels(srcR, srcG, srcB, dstR, dstG, dstB);

                if(allChannelFlags || channelFlags.testBit(red_pos))
                    dst[red_pos] = div(blend(src[red_pos], srcAlpha, dst[red_pos], dstAlpha, scale<channels_type>(dstR)), newDstAlpha);

                if(allChannelFlags || channelFlags.testBit(green_pos))
                    dst[green_pos] = div(blend(src[green_pos], srcAlpha, dst[green_pos], dstAlpha, scale<channels_type>(dstG)), newDstAlpha);

                if(allChannelFlags || channelFlags.testBit(blue_pos))
                    dst[blue_pos] = div(blend(src[blue_pos], srcAlpha, dst[blue_pos], dstAlpha, scale<channels_type>(dstB)), newDstAlpha);
            }

            return newDstAlpha;
        }
    }
};


template<class Traits, void compositeFunc(float, float, float, float&, float&, float&)>
class KoCompositeOpGenericHSL : public KoCompositeOpGenericHSLFunctor<Traits, detail::CompositeFunctionWrapperHSL<Traits, compositeFunc>>
{
protected:
    using base_class = KoCompositeOpGenericHSLFunctor<Traits, detail::CompositeFunctionWrapperHSL<Traits, compositeFunc>>;
public:
    using base_class::base_class;
};

/**
 * Generic CompositeOp for separable channel + alpha compositing functions
 *
 * A template to generate a KoCompositeOp class by just specifying a
 * blending/compositing function. This template works with compositing functions
 * for separable channels (means each channel of a pixel can be processed separately)
 * with taking alpha into consideration.
 * Note that because of special treating of alpha, any composite op function
 * needs to make alpha blending itself - the value of color that is written onto the projection
 * is the same that the composite function gives (compare with KoCompositeOpGenericHSL and KoCompositeOpGenericSC).
 */
template<class Traits, void compositeFunc(float, float, float&, float&), typename BlendingPolicy>
class KoCompositeOpGenericSCAlpha: public KoCompositeOpBase< Traits, KoCompositeOpGenericSCAlpha<Traits,compositeFunc,BlendingPolicy> >
{
    typedef KoCompositeOpBase< Traits, KoCompositeOpGenericSCAlpha<Traits,compositeFunc,BlendingPolicy> > base_class;
    typedef typename Traits::channels_type                                             channels_type;

    static const qint32 channels_nb = Traits::channels_nb;
    static const qint32 alpha_pos  = Traits::alpha_pos;

public:
    KoCompositeOpGenericSCAlpha(const KoColorSpace* cs, const QString& id, const QString& category)
        : base_class(cs, id, category) { }

public:
    template<bool alphaLocked, bool allChannelFlags>
    inline static channels_type composeColorChannels(const channels_type* src, channels_type srcAlpha,
                                                     channels_type*       dst, channels_type dstAlpha, channels_type maskAlpha,
                                                     channels_type opacity, const QBitArray& channelFlags)
    {
        using namespace Arithmetic;

        srcAlpha = mul(srcAlpha, maskAlpha, opacity);

        if(alphaLocked) {
            channels_type oldAlpha = dstAlpha;
            if(dstAlpha != zeroValue<channels_type>()) {
                for(qint32 i=0; i <channels_nb; i++) {
                    if(i != alpha_pos && (allChannelFlags || channelFlags.testBit(i))) {
                        float dstValueFloat = scale<float>(BlendingPolicy::toAdditiveSpace(dst[i]));
                        float dstAlphaFloat = scale<float>(oldAlpha);
                        compositeFunc(scale<float>(BlendingPolicy::toAdditiveSpace(src[i])), scale<float>(srcAlpha), dstValueFloat, dstAlphaFloat);
                        dst[i] = BlendingPolicy::fromAdditiveSpace(scale<channels_type>(dstValueFloat));
                    }
                }
            }

            return dstAlpha;
        }
        else {
            channels_type oldAlpha = dstAlpha;
            channels_type newDstAlpha = unionShapeOpacity(srcAlpha, dstAlpha);

            if(newDstAlpha != zeroValue<channels_type>()) {
                for(qint32 i=0; i <channels_nb; i++) {
                    if(i != alpha_pos && (allChannelFlags || channelFlags.testBit(i))) {
                        float dstFloat = scale<float>(BlendingPolicy::toAdditiveSpace(dst[i]));
                        float dstAlphaFloat = scale<float>(oldAlpha);
                        compositeFunc(scale<float>(BlendingPolicy::toAdditiveSpace(src[i])), scale<float>(srcAlpha), dstFloat, dstAlphaFloat);
                        dst[i] = BlendingPolicy::fromAdditiveSpace(scale<channels_type>(dstFloat));
                    }
                }
            }

            return newDstAlpha;
        }
    }
};

template<class Traits, void compositeFunc(float, float, float, float, float&, float&, float&)>
class KoCompositeOpGenericOVER: public KoCompositeOpBase< Traits, KoCompositeOpGenericOVER<Traits,compositeFunc> >
{
    typedef KoCompositeOpBase< Traits, KoCompositeOpGenericOVER<Traits,compositeFunc> > base_class;
    typedef typename Traits::channels_type                                             channels_type;

    static const qint32 red_pos   = Traits::red_pos;
    static const qint32 green_pos = Traits::green_pos;
    static const qint32 blue_pos  = Traits::blue_pos;

public:
    KoCompositeOpGenericOVER(const KoColorSpace* cs, const QString& id, const QString& category)
        : base_class(cs, id, category) { }

public:
    template<bool alphaLocked, bool allChannelFlags>
    inline static channels_type composeColorChannels(const channels_type* src, channels_type srcAlpha,
                                                     channels_type*       dst, channels_type dstAlpha, channels_type maskAlpha,
                                                     channels_type opacity, const QBitArray& channelFlags) {
        using namespace Arithmetic;
        srcAlpha = mul(srcAlpha, maskAlpha, opacity);

        if(srcAlpha == zeroValue<channels_type>()) return dstAlpha;

        channels_type newDstAlpha = (alphaLocked) ? dstAlpha : unionShapeOpacity(srcAlpha, dstAlpha);

        if(dstAlpha == zeroValue<channels_type>()) {
            if(allChannelFlags || channelFlags.testBit(red_pos))
                dst[red_pos] = src[red_pos];

            if(allChannelFlags || channelFlags.testBit(green_pos))
                dst[green_pos] = src[green_pos];

            if(allChannelFlags || channelFlags.testBit(blue_pos))
                dst[blue_pos] = src[blue_pos];
        } else {
            const float factor = float(srcAlpha) / newDstAlpha;

            float srcR = scale<float>(src[red_pos]);
            float srcG = scale<float>(src[green_pos]);
            float srcB = scale<float>(src[blue_pos]);

            float dstR = scale<float>(dst[red_pos]);
            float dstG = scale<float>(dst[green_pos]);
            float dstB = scale<float>(dst[blue_pos]);

            compositeFunc(srcR, srcG, srcB, 1.0f - factor, dstR, dstG, dstB);

            if(allChannelFlags || channelFlags.testBit(red_pos))
                dst[red_pos] = scale<channels_type>(dstR);

            if(allChannelFlags || channelFlags.testBit(green_pos))
                dst[green_pos] = scale<channels_type>(dstG);

            if(allChannelFlags || channelFlags.testBit(blue_pos))
                dst[blue_pos] = scale<channels_type>(dstB);
        }

        return newDstAlpha;
    }
};






template<class Traits, void compositeFunc(float, float, float, float, float&, float&, float&)>
class KoCompositeOpGenericCOPY: public KoCompositeOpBase< Traits, KoCompositeOpGenericCOPY<Traits,compositeFunc> >
{
    typedef KoCompositeOpBase< Traits, KoCompositeOpGenericCOPY<Traits,compositeFunc> > base_class;
    typedef typename Traits::channels_type                                             channels_type;

    static const qint32 red_pos   = Traits::red_pos;
    static const qint32 green_pos = Traits::green_pos;
    static const qint32 blue_pos  = Traits::blue_pos;

public:
    KoCompositeOpGenericCOPY(const KoColorSpace* cs, const QString& id, const QString& category)
        : base_class(cs, id, category) { }

public:
    template<bool alphaLocked, bool allChannelFlags>
    inline static channels_type composeColorChannels(const channels_type* src, channels_type srcAlpha,
                                                     channels_type*       dst, channels_type dstAlpha, channels_type maskAlpha,
                                                     channels_type opacity, const QBitArray& channelFlags) {
        using namespace Arithmetic;
        opacity = mul(maskAlpha, opacity);

        if (opacity == zeroValue<channels_type>()) return dstAlpha;

        channels_type newDstAlpha = (alphaLocked) ? dstAlpha : lerp(dstAlpha, srcAlpha, opacity);

        if (srcAlpha == zeroValue<channels_type>()) return newDstAlpha;

        if (dstAlpha == zeroValue<channels_type>()) {
            if(allChannelFlags || channelFlags.testBit(red_pos))
                dst[red_pos] = src[red_pos];

            if(allChannelFlags || channelFlags.testBit(green_pos))
                dst[green_pos] = src[green_pos];

            if(allChannelFlags || channelFlags.testBit(blue_pos))
                dst[blue_pos] = src[blue_pos];
        } else {
            /**
             * TODO: think about premultiplying the color channels before
             *       starting the blend (most probably, we don't need that,
             *       but who knows)
             */

            float srcR = scale<float>(src[red_pos]);
            float srcG = scale<float>(src[green_pos]);
            float srcB = scale<float>(src[blue_pos]);

            float dstR = scale<float>(dst[red_pos]);
            float dstG = scale<float>(dst[green_pos]);
            float dstB = scale<float>(dst[blue_pos]);

            compositeFunc(srcR, srcG, srcB, 1.0f - opacity, dstR, dstG, dstB);

            if(allChannelFlags || channelFlags.testBit(red_pos))
                dst[red_pos] = scale<channels_type>(dstR);

            if(allChannelFlags || channelFlags.testBit(green_pos))
                dst[green_pos] = scale<channels_type>(dstG);

            if(allChannelFlags || channelFlags.testBit(blue_pos))
                dst[blue_pos] = scale<channels_type>(dstB);
        }

        return newDstAlpha;
    }
};



#endif // _KOCOMPOSITEOP_GENERIC_H_
