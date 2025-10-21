#import <Foundation/Foundation.h>
#import <Metal/Metal.h>

@interface ProjectorClass : NSObject
   @property (nonatomic, assign) NSInteger no_norm;
   @property (nonatomic, strong) id<MTLDevice> device;
   @property (nonatomic, strong) NSMutableArray<id<MTLBuffer>> *d_Summ;
   @property (atomic, strong) id<MTLBuffer> d_output; // 2D/3D output
   @property (nonatomic, strong) id<MTLBuffer> d_im; // 3D volume
   @property (nonatomic, strong) NSMutableArray<id<MTLBuffer>> *d_rhs_os;
   //@property (nonatomic, strong) NSMutableArray<id<MTLBuffer>> *d_norm; // Normalization data
   //@property (nonatomic, strong) NSMutableArray<id<MTLBuffer>> *d_scat; // Scatter data
   @property (nonatomic, strong) id<MTLBuffer> d_attenB; // Attenuation image
   @property (nonatomic, strong) id<MTLBuffer> d_maskBP; // 2D/3D BP mask

- (NSInteger)addProjector:(scalarStruct)inputScalars
                           weighting:(Weighting)w_vec
                              //method:(RecMethodsBox *)methodList
                    headerDirectory:(NSString *)headerDirectory
                                type:(NSInteger)type;
                              
- (NSInteger)createBuffers:(scalarStruct)inputScalars
                           weighting:(Weighting)w_vec
                           x:(const float *)x
                       zDet:(const float *)z_det
                   xyIndex:(const uint32_t *)xy_index
                    zIndex:(const uint16_t *)z_index
                    L:(const uint16_t *)L
                    pituus:(const int64_t *)pituus
                     atten:(const float *)atten
                      norm:(const float *)norm
                 extraCorr:(const float *)extraCorr
                    length:(const int64_t *)length
                    type:(NSUInteger)type;
@end