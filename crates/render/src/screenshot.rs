// ---------------------------------------------------------------------------
// Screenshot — capture offscreen framebuffer to PNG
// ---------------------------------------------------------------------------

use std::io::Write;
use std::path::Path;

use crate::field_pipeline::FieldPipeline;

/// Capture the current offscreen framebuffer as PNG bytes.
///
/// This blocks until the GPU readback is complete. Use only for
/// user-triggered screenshot actions.
pub fn capture_screenshot(
    device: &wgpu::Device,
    queue: &wgpu::Queue,
    pipeline: &FieldPipeline,
) -> Option<Vec<u8>> {
    let [width, height] = pipeline.framebuffer_size();
    if width == 0 || height == 0 {
        return None;
    }

    let bytes_per_row = (width * 4 + 255) & !255; // align to 256

    // Create staging buffer for readback
    let staging = device.create_buffer(&wgpu::BufferDescriptor {
        label: Some("screenshot-staging"),
        size: (bytes_per_row * height) as u64,
        usage: wgpu::BufferUsages::MAP_READ | wgpu::BufferUsages::COPY_DST,
        mapped_at_creation: false,
    });

    // Copy texture to buffer
    let mut encoder = device.create_command_encoder(&wgpu::CommandEncoderDescriptor {
        label: Some("screenshot-encoder"),
    });
    encoder.copy_texture_to_buffer(
        wgpu::TexelCopyTextureInfo {
            texture: pipeline.color_texture(),
            mip_level: 0,
            origin: wgpu::Origin3d::ZERO,
            aspect: wgpu::TextureAspect::All,
        },
        wgpu::TexelCopyBufferInfo {
            buffer: &staging,
            layout: wgpu::TexelCopyBufferLayout {
                offset: 0,
                bytes_per_row: Some(bytes_per_row),
                rows_per_image: Some(height),
            },
        },
        wgpu::Extent3d {
            width,
            height,
            depth_or_array_layers: 1,
        },
    );
    queue.submit(Some(encoder.finish()));

    // Map the buffer and wait
    let slice = staging.slice(..);
    let (tx, rx) = std::sync::mpsc::channel();
    slice.map_async(wgpu::MapMode::Read, move |result| {
        let _ = tx.send(result);
    });
    device.poll(wgpu::PollType::wait_indefinitely()).unwrap();

    if rx.recv().ok()?.is_err() {
        return None;
    }

    // Read pixels and encode PNG
    let data = slice.get_mapped_range();
    let mut png_data = Vec::new();
    encode_png(&data, width, height, bytes_per_row, &mut png_data).ok()?;
    drop(data);
    staging.unmap();

    Some(png_data)
}

/// Save screenshot to a file path.
pub fn save_screenshot(path: &Path, png_bytes: &[u8]) -> std::io::Result<()> {
    std::fs::write(path, png_bytes)
}

/// Encode RGBA8 pixel data as PNG.
fn encode_png(
    data: &[u8],
    width: u32,
    height: u32,
    src_bytes_per_row: u32,
    output: &mut Vec<u8>,
) -> std::io::Result<()> {
    // Minimal PNG encoder without external crate dependency
    // PNG signature
    output.write_all(&[137, 80, 78, 71, 13, 10, 26, 10])?;

    // IHDR chunk
    let mut ihdr = Vec::new();
    ihdr.extend_from_slice(&width.to_be_bytes());
    ihdr.extend_from_slice(&height.to_be_bytes());
    ihdr.push(8); // bit depth
    ihdr.push(6); // color type: RGBA
    ihdr.push(0); // compression
    ihdr.push(0); // filter
    ihdr.push(0); // interlace
    write_png_chunk(output, b"IHDR", &ihdr)?;

    // IDAT chunk: build raw pixel data with filter bytes
    let row_bytes = (width * 4) as usize;
    let mut raw = Vec::with_capacity((1 + row_bytes) * height as usize);
    for y in 0..height as usize {
        raw.push(0); // no filter
        let src_offset = y * src_bytes_per_row as usize;
        raw.extend_from_slice(&data[src_offset..src_offset + row_bytes]);
    }

    // Compress with deflate (use miniz_oxide via flate2 if available, else store)
    let compressed = deflate_compress(&raw);
    write_png_chunk(output, b"IDAT", &compressed)?;

    // IEND chunk
    write_png_chunk(output, b"IEND", &[])?;

    Ok(())
}

fn write_png_chunk(output: &mut Vec<u8>, chunk_type: &[u8; 4], data: &[u8]) -> std::io::Result<()> {
    output.write_all(&(data.len() as u32).to_be_bytes())?;
    output.write_all(chunk_type)?;
    output.write_all(data)?;
    // CRC32 of type + data
    let crc = crc32(chunk_type, data);
    output.write_all(&crc.to_be_bytes())?;
    Ok(())
}

fn crc32(chunk_type: &[u8], data: &[u8]) -> u32 {
    let mut crc: u32 = 0xFFFF_FFFF;
    for &byte in chunk_type.iter().chain(data.iter()) {
        crc ^= byte as u32;
        for _ in 0..8 {
            if crc & 1 != 0 {
                crc = (crc >> 1) ^ 0xEDB8_8320;
            } else {
                crc >>= 1;
            }
        }
    }
    crc ^ 0xFFFF_FFFF
}

/// Simple deflate compression using zlib-compatible format (stored blocks).
/// This produces a valid but uncompressed zlib stream.
fn deflate_compress(data: &[u8]) -> Vec<u8> {
    let mut output = Vec::new();
    // Zlib header: CM=8 (deflate), CINFO=7, FCHECK
    output.push(0x78);
    output.push(0x01);

    // Split data into stored blocks (max 65535 bytes each)
    let mut offset = 0;
    while offset < data.len() {
        let remaining = data.len() - offset;
        let block_size = remaining.min(65535);
        let is_final = offset + block_size >= data.len();

        output.push(if is_final { 0x01 } else { 0x00 });
        let len = block_size as u16;
        let nlen = !len;
        output.extend_from_slice(&len.to_le_bytes());
        output.extend_from_slice(&nlen.to_le_bytes());
        output.extend_from_slice(&data[offset..offset + block_size]);
        offset += block_size;
    }

    // Adler32 checksum
    let adler = adler32(data);
    output.extend_from_slice(&adler.to_be_bytes());
    output
}

fn adler32(data: &[u8]) -> u32 {
    let mut a: u32 = 1;
    let mut b: u32 = 0;
    for &byte in data {
        a = (a + byte as u32) % 65521;
        b = (b + a) % 65521;
    }
    (b << 16) | a
}
