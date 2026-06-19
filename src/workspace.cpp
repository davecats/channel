#include "channel/workspace.hpp"

#include "channel/memory.hpp"

#include <algorithm>

namespace channel {

WorkspaceLease::WorkspaceLease(WorkspaceArena* arena, std::string owner, std::size_t bytes, std::uint8_t* base)
    : arena_(arena), owner_(std::move(owner)), bytes_(bytes), base_(base) {}

WorkspaceLease::WorkspaceLease(WorkspaceLease&& other) noexcept
    : arena_(other.arena_),
      owner_(std::move(other.owner_)),
      bytes_(other.bytes_),
      cursor_(other.cursor_),
      base_(other.base_) {
  other.arena_ = nullptr;
  other.bytes_ = 0;
  other.cursor_ = 0;
  other.base_ = nullptr;
}

WorkspaceLease& WorkspaceLease::operator=(WorkspaceLease&& other) noexcept {
  if (this != &other) {
    release();
    arena_ = other.arena_;
    owner_ = std::move(other.owner_);
    bytes_ = other.bytes_;
    cursor_ = other.cursor_;
    base_ = other.base_;
    other.arena_ = nullptr;
    other.bytes_ = 0;
    other.cursor_ = 0;
    other.base_ = nullptr;
  }
  return *this;
}

WorkspaceLease::~WorkspaceLease() {
  release();
}

std::size_t WorkspaceLease::align_offset(std::size_t offset, std::size_t alignment) {
  const std::size_t rem = offset % alignment;
  return rem == 0 ? offset : offset + alignment - rem;
}

void WorkspaceLease::release() {
  if (arena_) {
    arena_->release(owner_);
    arena_ = nullptr;
  }
}

WorkspaceArena::WorkspaceArena(WorkspaceBackend backend) : backend_(backend) {
  if (backend_ == WorkspaceBackend::UmpireSpace) {
    initialize_memory_pool();
  }
}

WorkspaceArena::~WorkspaceArena() = default;

WorkspaceLease WorkspaceArena::lease(std::size_t bytes, const std::string& owner) {
  if (owned_) {
    throw std::runtime_error("WorkspaceArena already owned by " + owner_);
  }
  reserve(bytes);
  owned_ = true;
  owner_ = owner;
  high_water_ = std::max(high_water_, bytes);
  return WorkspaceLease(this, owner, bytes, storage_data());
}

void WorkspaceArena::reserve(std::size_t bytes) {
  if (owned_) {
    throw std::runtime_error("WorkspaceArena cannot reserve while owned by " + owner_);
  }
  if (bytes <= capacity_) return;

  const std::size_t alloc_bytes = bytes + alignment;
  storage_.resize("workspace_arena", alloc_bytes);
  capacity_ = bytes;
  high_water_ = std::max(high_water_, bytes);
}

void WorkspaceArena::release(const std::string& owner) {
  if (!owned_) {
    throw std::runtime_error("WorkspaceArena release requested while unowned");
  }
  if (owner != owner_) {
    throw std::runtime_error("WorkspaceArena owned by " + owner_ + ", release requested by " + owner);
  }
  owned_ = false;
  owner_.clear();
}

std::uint8_t* WorkspaceArena::storage_data() {
  return storage_.data();
}

} // namespace channel
